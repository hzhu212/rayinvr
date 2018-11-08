Some part of rayinvr is too old to be compiled correctly. This project lists some patches to make rayinvr run properly on Ubuntu Linux 18.04.

# How to compile

Just follow the `Makefile` under each subdirectory. But before that, you should make some patches first.

## Install tool chains

```sh
apt install gfortran libx11-dev
```

## Patches for compiling on Ubuntu 18.04

### Makefile: unilink not found

We can just let `unilink` go, and use `$(FC)` instead. Open the `Makefile` that throws this error, and replace the line:

```sh
unilink ${TRAMP_OBJS}
```

with

```sh
$(FC) -o main ${TRAMP_OBJS}
```

> `$(FC)` is an implicit variable of `make`, which sets the command of Fortran Compliler. It equals 'f77' by default, but you can set it temperarily when running `make` command. For example: `make ... FC=gfortran`.
> Similarly, you can set variable `CC` to specify your C compiler.

### Issue of fortran "Variable FORMAT expression"

`gfortran` does not support "Variable FORMAT expression". For example, `pltlib.f:198` reads:


```fortran
5     format(i2/i10/4e15.5,<nchar>a1)
```

The expression `<nchar>` is a "Variable FORMAT expression". We can find out from the code that `nchar` is an local integer variable. When compiling this source file, `gfortran` will throw a compiling error:

```txt
gfortran Error: Unexpected element ‘<’ in format string at (1)
```

To fix this issue, there are 2 common ways:

1. As GUN Fortran recommended, we can construct a format string dynamicly. See [here](http://gcc.gnu.org/onlinedocs/gfortran/Variable-FORMAT-expressions.html)

2. Simply replace `<nchar>` with a number which is large enough, say `10000`. Now, the `format` line will looks like this:

    ```fortran
    5     format(i2/i10/4e15.5,10000a1)
    ```

### Issue of character '$'

Using '$' in variable names is not allowed by default in gfortran. Or it will throw the following error:

```txt
# when building xrayinvr
# in rayinvr/trc.f:498:17:
Fatal Error: Invalid character ‘$’ at (1). Use ‘-fdollar-ok’ to allow it as an extension compilation terminated
```

To fix this error, you should add '-fdollar-ok' option to `f77` command. However, we don't use `f77` command directly; what we need to do is change the `Makefile` a bit:

```make
FFLAGS = -O -fdollar-ok
```

### Change backgroud color of Xbuplot

`pltlib/xbuplot.c` used `Xlib` to create window, the foreground and background settings are in line 99, 100:

```c
  /*  create opaque window  */
  win = XCreateSimpleWindow (display, RootWindow (display, screen_num),
                x, y, win_width, win_height, border_width,
        WhitePixel (display, screen_num),
        BlackPixel (display, screen_num));
```

### Read v.in with higher precision

The default precision of `v.in` is `f7.2`. If you want to read a higher precision version of `v.in`(say, `f8.3`), just open the `main.f` that reads the `v.in`, find and change the lines like

```fortran
15       format(i2,1x,10f7.2)
25       format(3x,10i7)
25       format(3x,10(5x,i2))
```

into

```fortran
15       format(i2,1x,10f8.3)
25       format(3x,10i8)
25       format(3x,10(6x,i2))
```

### Some other syntax errors

`pltsyn/pltsec.f`:253, change from:

```fortran
108     if(namp.eq.0) write(*,*, fmt="(/
```

into:

```fortran
108     if(namp.eq.0) write(*, fmt="(/
```



## An example for compiling: `tramp`

Assuming we are in rayinvr source directory now. Firstly, go to `tramp` subdirectory:

```sh
cd tramp
```

Secondly, modify the `Makefile` and set proper input and output directory. ENSURE that the output directory exists. Then you can start to make:

```sh
make tramp
make xtramp
```
