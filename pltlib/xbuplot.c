#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <stdio.h>

#include "xbuplot_icon"

#define PIXELS_PER_INCH (100)
#define MM_PER_INCH (25.4)
#define XFLUSH 0
#define MAX_COLOURS 16

/*  These are used as arguments to nearly every Xlib routine, so it saves 
 *  routine arguments to declare them global.  If there were 
 *  additional source files, they would be declared extern there. 
 */
Display *display;
int screen_num;
Window win;
GC gc;
XFontStruct *font_info;
unsigned int win_width, win_height;  /*  window size  */
float win_width_mm, win_height_mm;  /*  window size in mm  */
int last_x = 0;  /*  last x position of pen  */
int last_y = 0;  /*  last y position of pen  */
static int can_use_colours;
static unsigned long colour_value [MAX_COLOURS];
static char *colour_name[] = {"white", "black", "red", "green",
			"blue", "cyan", "magenta", "yellow",
			"orange", "yellow green", "turquoise", "blue violet",
			"violet", "purple", "lavender", "aquamarine"};

void xbuinit_ ();
void xbuplot_ ();
void xbuline_ ();
void xbulines_ ();
void xbusymbol_ ();
void xbuerase_ ();
void xbuflush_ ();
void xbuclose_ ();
unsigned int mm_to_pixels ();
float pwr ();
int colour_available ();
unsigned long get_colour ();
void xbucolour_ ();
void xburect_ ();
void xbuevent_ ();
void xbucircle_ ();


/********  xbuinit  ********/
void xbuinit_ (x_mm_size, y_mm_size)
  float *x_mm_size, *y_mm_size;
{
  unsigned int x, y;  /*  window position  */
  unsigned int border_width = 4;  /*  four pixels  */
  XEvent event;
  char *display_name = NULL;
  char *font_name = "6x12";
  Pixmap icon_pixmap;
  char *icon_name = "Xbuplot";
  char *window_name = "Xbuplot";
  XTextProperty windowName, iconName;
  XSizeHints size_hints;
  XWMHints wm_hints;
  XClassHint class_hints;
  XSetWindowAttributes set_window_attributes;
  int exposed_once = 0;  /*  flag if window has been exposed yet  */
  
  /*  connect to X server  */
  if ((display = XOpenDisplay (display_name)) == NULL)
  {
    (void) fprintf (stderr, "ERROR: cannot connect to X server %s\n",
			XDisplayName (display_name));
    exit (-1);
  }

  /*  get screen_num from DefaultScreen  */
  screen_num = DefaultScreen (display);

  /*  Note that in a real application, x and y would default to 0
   *  but would be settable from the command line or resource database.  
   */
  x = y = 0;

  /*  save the window size in mm  */
  win_width_mm = *x_mm_size;
  win_height_mm = *y_mm_size;

  /*  convert window size from mm to pixels  */
  win_width = mm_to_pixels (*x_mm_size);
  win_height = mm_to_pixels (*y_mm_size);

  /*  create opaque window  */
  win = XCreateSimpleWindow (display, RootWindow (display, screen_num), 
				x, y, win_width, win_height, border_width,
				WhitePixel (display, screen_num), 
				BlackPixel (display, screen_num));

  /*  Create pixmap of depth 1 (bitmap) for icon  */
  icon_pixmap = XCreateBitmapFromData (display, win, xbuplot_icon_bits, 
					xbuplot_icon_width,xbuplot_icon_height);

  /*  x, y, width, and height hints are now taken from
   *  the actual settings of the window when mapped. Note
   *  that PPosition and USSize must be specified anyway. 
   *  We specify PMaxSize and PMinSize so the user can't resize the window
   */
  size_hints.flags = PPosition | USSize | PMaxSize | PMinSize;
  size_hints.min_width = size_hints.max_width = win_width;
  size_hints.min_height = size_hints.max_height = win_height;

  /*  These calls store window_name and icon_name into XTextProperty 
   *  structures and set their other fields properly. 
   */
  if (XStringListToTextProperty (&window_name, 1, &windowName) == 0)
  {
    (void) fprintf (stderr, 
			"ERROR: structure allocation for windowName failed.\n");
    exit(-1);
  }
  if (XStringListToTextProperty (&icon_name, 1, &iconName) == 0)
  {
    (void) fprintf (stderr, 
			"ERROR: structure allocation for iconName failed.\n");
    exit(-1);
  }

  /*  set window manager hints  */
  wm_hints.initial_state = NormalState;
  wm_hints.input = True;
  wm_hints.icon_pixmap = icon_pixmap;
  wm_hints.flags = StateHint | IconPixmapHint | InputHint;

  class_hints.res_name = "xbuplot";
  class_hints.res_class = "xbuplot";
  
  XSetWMProperties (display, win, &windowName, &iconName, NULL, 0, 
			&size_hints, &wm_hints, &class_hints);

  /*  Set backing store to on since we can't service Expose events  */
  set_window_attributes.backing_store = Always;
  XChangeWindowAttributes (display,win, CWBackingStore, &set_window_attributes);

  /*  Select event types wanted  */
  XSelectInput (display, win, ExposureMask | StructureNotifyMask | 
		ButtonPressMask | KeyPressMask);

  /*  Load font and get font information structure  */
  if ((font_info = XLoadQueryFont (display, font_name)) == NULL)
  {
    (void) fprintf (stderr, "ERROR: Cannot open %s font\n", font_name);
    exit( -1 );
  }

  /*  Create Graphics Context and use defaults XGCvalues  */
  gc = XCreateGC (display, win, 0, NULL);

  /*  specify font  */
  XSetFont (display, gc, font_info->fid);

  /*  check to see if colours are available  */
  /*  if not set foreground colour to white  */
  if (!(can_use_colours = colour_available ()))
    XSetForeground (display, gc, WhitePixel (display, screen_num));

  /*  Display window  */
  XMapWindow (display, win);

}


/********  xbuplot_  ********/
void xbuplot_ (x_mm, y_mm, pen_down)
  float *x_mm, *y_mm;  /*  position in mm of where to move pen  */
  int *pen_down;  /*  pen = 1 --> pen down; pen = 0 --> pen up  */
{
  int x, y;  /*  position in pixels of where to move pen  */

  x = mm_to_pixels (*x_mm);
  y = win_height - mm_to_pixels (*y_mm);
  if (*pen_down)
  {
    XDrawLine (display, win, gc, last_x, last_y, x, y);
#if XFLUSH
    XFlush (display);
#endif
  }
  last_x = x;
  last_y = y;
}
  

/********  xbuline_  ********/
void xbuline_ (x1_mm, y1_mm, x2_mm, y2_mm)
  float *x1_mm, *y1_mm;  /*  first point in line  */
  float *x2_mm, *y2_mm;  /*  second point in line  */
{
  int x1, y1, x2, y2;  /*  x y coords of line in pixels  */

  x1 = mm_to_pixels (*x1_mm);
  y1 = win_height - mm_to_pixels (*y1_mm);
  x2 = mm_to_pixels (*x2_mm);
  y2 = win_height - mm_to_pixels (*y2_mm);
  XDrawLine (display, win, gc, x1, y1, x2, y2);
#if XFLUSH
  XFlush (display);
#endif
  last_x = x2;
  last_y = y2;
}
  

/********  xbulines_  ********/
void xbulines_ (x_mm, y_mm, num_points)
  float *x_mm, *y_mm;  /*  points in line  */
  int *num_points;  /*  number of points in line  */
{
  int x, y;  /*  x y coords of line in pixels  */
  int point = 0;

  last_x = mm_to_pixels (x_mm [0]);
  last_y = win_height - mm_to_pixels (y_mm [0]);
  for (point = 1; point <= *num_points; point++)
  {
    x = mm_to_pixels (x_mm [point - 1]);
    y = win_height - mm_to_pixels (y_mm [point - 1]);
    XDrawLine (display, win, gc, last_x, last_y, x, y);
    last_x = x;
    last_y = y;
  }
#if XFLUSH
  XFlush (display);
#endif
}

  
/********  xbusymbol_  ********/
void xbusymbol_ (x_mm, y_mm, height_mm, string, angle, num_chars)
  float *x_mm, *y_mm;  /*  point to draw string at  */
  float *height_mm;  /*  height of text (not yet implemented)  */
  char *string;  /*  string to be drawn  */
  float *angle;  /*  angle in degrees of text. only 0 and 90 supported  */
  int *num_chars;  /*  number of characters in string  */
{
  int x, y;  /*  x and y coords of the bottom left corner of text  */
  int i;
  int height;
  char character [2];

  x = mm_to_pixels (*x_mm);
  y = win_height - mm_to_pixels (*y_mm);
  height = mm_to_pixels (*height_mm);
  character [1] = 0;

  if ((int) (*angle + 0.5) == 0) 
  {
    XDrawString (display, win, gc, x, y, string, *num_chars);
  } 
  else if ((int) (*angle + 0.5) == 90)
  {
    for (i = 1; i <= *num_chars; i++)
    {
      character [0] = string [i - 1];
      XDrawString (display, win, gc, x, 
			(unsigned int) (y - ((*num_chars - i) * 1.2 * height)), 
			character, 1);
    }
  }

#if XFLUSH
  XFlush (display);
#endif
}


/********  xbunumber_  ********/
void xbunumber_ (x_mm, y_mm, height_mm, number, angle, num_dec_points)
  float *x_mm, *y_mm;  /*  point to draw string at  */
  float *height_mm;  /*  height of text (not yet implemented)  */
  float *number;  /*  number to be plotted  */
  float *angle;   /*  angle in degrees of text, only 0 and 90 supported  */
  int *num_dec_points;  /*  number of digits to follow decimal point  */
{
  int x, y;  /*  x and y coords of the bottom left corner of text  */
  char string[80];
  int num_chars;
  float rounded_number;
  int ndp;
  int i;
  int height;
  char character [2];

  x = mm_to_pixels (*x_mm);
  y = win_height - mm_to_pixels (*y_mm);
  height = mm_to_pixels (*height_mm);
  character [1] = 0;
  ndp = *num_dec_points;
  if (*number > 0.0)
  {
    rounded_number = (float) ((long) (*number * pwr (10.0, ndp) + 0.5) 
                        / pwr (10.0,  ndp));
  }
  else
  {
    rounded_number = (float) ((long) (*number * pwr (10.0, ndp) - 0.5) 
                        / pwr (10.0,  ndp));
  }
  sprintf (string, "%g", rounded_number);
  num_chars = strlen (string);

  if ((int) (*angle + 0.5) == 0) 
  {
    XDrawString (display, win, gc, x, y, string, num_chars);
  } 
  else if ((int) (*angle + 0.5) == 90)
  {
    for (i = 1; i <= num_chars; i++)
    {
      character [0] = string [i - 1];
      XDrawString (display, win, gc, x, 
                        (unsigned int) (y - ((num_chars - i) * 1.2 * height)), 
                        character, 1);
    }
  }

#if XFLUSH
  XFlush (display);
#endif
}


/********  xbuerase_  ********/
void xbuerase_ ()
{

  XClearArea (display, win, 0, 0, 0, 0, False);
#if XFLUSH
  XFlush (display);
#endif
}


/********  xbuflush_  ********/
void xbuflush_ ()
{

  XFlush (display);
}


/********  xbuclose_  ********/
void xbuclose_ ()
{
        XUnloadFont (display, font_info->fid);
        XFreeGC (display, gc);
        XCloseDisplay (display);
}


/********  mm_to_pixels  ********/
unsigned int mm_to_pixels (mm)
  float mm;
{
  return ((unsigned int) (mm * PIXELS_PER_INCH / MM_PER_INCH + 0.5));
}


/********  pwr  ********/
float pwr (number, exponent)
  float number;
  int exponent;
{
  int i;
  float value;

  value = 1.0;
  
  for (i = 1; i <= exponent; i++)
  {
    value *= number;
  }
  return (value);
}


/********  colour_available  ********/
/*  find out in the display supports colour. 
 *  i.e. PseudoColor, TrueColor, DirectColor, and StaticColor
 *  return 0 is colour_available == false
 *  return 1 is colour_available == true
 */
int colour_available ()
{
  int default_depth;
  int class = 5;
  XVisualInfo visual_info;

  default_depth = DefaultDepth (display, screen_num);
  if (default_depth == 1) 
  {
    /*  Must be static grey. use black and white  */
    /*  return value of 0 for "colour_available == false"  */
    return (0);
  }

  /*  check to see if a visual class that supports colour is available  */
  while (!XMatchVisualInfo (display, screen_num, default_depth, class--, 
				&visual_info))
    ;
  class++;  /*  set class back to which class was found  */
  if (class < 2)  /*  colour visual classes are 2 to 5  */
  {
    /*  No color visual available at default_depth. use black and white  */
    /*  return value of 0 for "colour_available == false"  */
    return (0);
  }
  /*  color visual is available. colours are available  */
  /*  return value of 1 for "colour_available == true"  */
  return (1);
}


/********  get_colour  ********/
/*  returns pixel value for given colour  */
/*  allocates colour if necessary  */
unsigned long get_colour (colour)
  int colour;
{
  Colormap default_colour_map;
  XColor exact_def;

  /*  check to make sure colour is in range, if not return white  */
  if (colour > MAX_COLOURS - 1)
    return (WhitePixel (display, screen_num));

  /*  check to see if colour has alread been allocated. if so return value  */
  if (colour_value [colour])  /*  colour_value is non-zero if allocated  */
    return (colour_value [colour]);

  /*  get default colour map  */
  default_colour_map   = DefaultColormap (display, screen_num);

  /*  parse colour name to see if colour is available. if not use white  */
  if (!XParseColor (display, default_colour_map, colour_name [colour], 
			&exact_def)) 
    return (WhitePixel (display, screen_num));

  /*  allocate colour. if not available use white  */
  if (!XAllocColor (display, default_colour_map, &exact_def))
    return (WhitePixel (display, screen_num));

  /*  we got the colour so store it and return it's pixel value  */
  colour_value [colour] = exact_def.pixel;
  return (exact_def.pixel);
}


/********  xbucolour  ********/
/*  set the foreground colour in the gc to the value of colour  */
void xbucolour_ (colour)
  int *colour;
{

  if (can_use_colours) 
    XSetForeground (display, gc, get_colour (*colour));
}


/********  xburect  ********/
/*  draw a rectangle and, possibly, fill it  */
void xburect_ (x1_mm, y1_mm, x2_mm, y2_mm, fill)
  float *x1_mm, *y1_mm;  /*  bottom left corner in mm  */
  float *x2_mm, *y2_mm;  /*  top right corner in mm  */
  int *fill;  /*  fill flag: 0 == no fill, 1 == fill  */
{
  unsigned int x1, y1;  /*  bottom left corner in pixels  */
  unsigned int x2, y2;  /*  top right right corner in pixels  */

  x1 = mm_to_pixels (*x1_mm);
/* MA following line changed from y1 = win_height - mm_to_pixels (*y1_mm) */
  y1 = mm_to_pixels (*y1_mm);
  x2 = mm_to_pixels (*x2_mm);
/* MA following line changed from y2 = win_height - mm_to_pixels (*y2_mm) */
  y2 = mm_to_pixels (*y2_mm);

  if (*fill)
  {
    XDrawRectangle (display, win, gc, x1, y1, x2 - x1, y2 - y1);
    XFillRectangle (display, win, gc, x1 + 1, y1 + 1, x2 - x1 - 1, y2 - y1 - 1);
  }
  else
  {
    XDrawRectangle (display, win, gc, x1, y1, x2 - x1, y2 - y1);
    XFillRectangle (display, win, gc, x1 + 1, y1 + 1, x2 - x1 - 1, y2 - y1 - 1);
  }
}


/********  xbuevent  ********/
/*  wait for a button or key press event  */
void xbuevent_ (x_mm, y_mm, button, key)
  float *x_mm, *y_mm;  /*  position of cursor in mm  */
  int *button;  /*  button that was pressed  */
  char *key;  /*  key that was pressed  */
{
  XEvent event;
  KeySym keysym;
  XComposeStatus compose;
  char buffer = 0;
  int buflen = 1;
  int got_event = 0;

  *button = *key = 0;	/* initialize button and key to zero */

  /* loop until we get a key or button press that we want */
  while (!got_event)
  {
    /* get next button or key press event */
    XMaskEvent (display, ButtonPressMask | KeyPressMask, &event);
    switch (event.type)
    {
      case ButtonPress:
        *button = event.xbutton.button;
        *x_mm = (double) event.xbutton.x / (double) win_width * win_width_mm;
        *y_mm = (double) (win_height - event.xbutton.y) / (double) win_height 
		* win_height_mm;
        got_event++;
        break;
      case KeyPress:
        XLookupString (&event, &buffer, buflen, &keysym, &compose);
        switch (keysym)
        {
          case XK_Shift_L:	/* ignore a key press from a key modifier */
          case XK_Shift_R:
          case XK_Control_L:
          case XK_Control_R:
          case XK_Caps_Lock:
          case XK_Shift_Lock:
          case XK_Meta_L:
          case XK_Meta_R:
          case XK_Alt_L:
          case XK_Alt_R:
          case XK_Super_L:
          case XK_Super_R:
          case XK_Hyper_L:
          case XK_Hyper_R:
            break;
          default:
            *key = buffer;
            *x_mm = (double) event.xkey.x / (double) win_width * win_width_mm;
            *y_mm = (double) (win_height - event.xkey.y) / (double) win_height 
			* win_height_mm;
            got_event++;
            break;
        }
    }
  }
}


/********  xbucircle  ********/
/*  draw a circle  */
void xbucircle_ (x_mm, y_mm, size_mm, fill_flag)
  float *x_mm, *y_mm;  /*  position of cursor in mm  */
  float *size_mm;  /*  size of circle  */
  int *fill_flag;  /*  fill flag: 0 = no fill, 1 = fill  */
{
  int x, y;  /*  x y coords of center of circle in pixels  */
  int size;  /*  size of circle in pixels  */

  x = mm_to_pixels (*x_mm);
  y = win_height - mm_to_pixels (*y_mm);
  size = mm_to_pixels (*size_mm);

  XDrawArc (display, win, gc, x, y, size, size, 0, 64 * 360);
  if (*fill_flag) XFillArc (display, win, gc, x, y, size, size, 0, 64 * 360);
#if XFLUSH
  XFlush (display);
#endif
}
