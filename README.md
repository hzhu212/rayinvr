A modified version of [rayinvr](http://terra.rice.edu/department/faculty/zelt/rayinvr.html) that can work properly on modern Linux(e.g. Ubuntu 18.04) and Windows(e.g. Windows 10) platform.

# How to use

## Prepare tool-chains

On Linux, you need to install `x11` library and `gfortran`. For example, on Ubuntu 18.04:

```sh
sudo apt install gfortran libx11-dev
```

On Windows, you need to install [GNU Make](http://gnuwin32.sourceforge.net/packages/make.htm) and [mingw-w64](http://www.mingw-w64.org/doku.php) manually.

## Build

Assuming we are in root directory of rayinvr source now. 

Inside each subdirectory there is a `Makefile` to make two versions of each program. For example, you can build `tramp` as the following steps:

```sh
cd tramp
make xtramp FC=gfortran CC=gcc
make tramp FC=gfortran CC=gcc
# or just `make FC=gfortran CC=gcc`, will build both the two versions.
```

PS: the keyword arguments `FC=gfortran CC=gcc` means setting Fortran Compiler to `gfortran` and C Compiler to `gcc`. By default, the two variables are `f77` and `cc` respectively.

The first version, prefixed by the letter `x`, is based on `X11` graphics system. The second version is originally for users who have access to the commercial graphics package called Uniras. However, the second version would also be helpful if you just need to calculate without ploting. On Windows, there is no `x11` library so that the second version is the only choice.

After this step, we will get `rayinvr.exe`, `tramp.exe` etc. under `build` directory.

## Add to PATH

**This step is optional.**

Add `build` directory to `PATH`, then we will be able to run `rayinvr` and other commands at any path.

## Run examples

In original rayinvr repository, Zelt provided 8 example datasets, which is helpful to test our build results. 

But we can't run our commands on these datasets directly, because I have changed `main.f` to handle higher float precision when reading `v.in` file (see [below](#read-vin-with-higher-precision)). If you have tried running commands at `examples` directory, you will probably end up with an error msg like this:

```
At line 236 of file main.f (unit = 20, file = 'v.in')
Fortran runtime error: Bad value during floating point read

Error termination. Backtrace:
```

But don't worry, I have made a python script which can convert old style v.in file into high precision version, and I have put all converted example datasets in `examples-new` directory.

Now, let's try example3 (Why not example1? because example3 is neither too simple, nor too complicated):

```cmd
# On Windows 10 WSL for my case

# go to dataset directory
cd examples-new/e3

# run rayinvr.exe
../../build/rayinvr.exe
```

The output will be as below, basically consistent with that of original rayinvr:

```
shot#   1:   ray code  2.2:    59 rays traced
shot#   1:   ray code  3.2:    60 rays traced
shot#  -2:   ray code  2.2:    17 rays traced
shot#  -2:   ray code  3.2:     7 rays traced
shot#   2:   ray code  2.2:    63 rays traced
shot#   2:   ray code  3.2:    73 rays traced
shot#  -3:   ray code  2.2:    34 rays traced
shot#  -3:   ray code  3.2:    26 rays traced
shot#   3:   ray code  2.2:    65 rays traced
shot#   3:   ray code  3.2:    73 rays traced

|----------------------------------------------------------------|
|                                                                |
| total of   2089 rays consisting of    21239 points were traced |
|                                                                |
|           model consists of  4 layers and  40 blocks           |
|                                                                |
|----------------------------------------------------------------|


Number of data points used:      477
RMS traveltime residual:       0.068
Normalized chi-squared:       46.331

Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG IEEE_DIVIDE_BY_ZERO IEEE_OVERFLOW_FLAG IEEE_DENORMAL
```

Now, start your work with rayinvr!

# Problems in original rayinvr and fix

This section lists what errors the original rayinvr throws and how I fix them.

## Makefile: unilink not found

Just abandon `unilink` command and use `$(FC)` instead. Open the `Makefile` that throws this error, and change the line:

```sh
unilink ${TRAMP_OBJS}
```

to

```sh
$(FC) -o main ${TRAMP_OBJS}
```

PS: `$(FC)` will be replaced with your Fortran Compiler, for example `gfortran`.

## Issue of fortran "Variable FORMAT expression"

Both `f77` and `gfortran` do not support "Variable FORMAT expression", which is used in several rayinvr programs. For example, `pltlib.f:198` reads:

```fortran
5     format(i2/i10/4e15.5,<nchar>a1)
```

This line includes a "Variable FORMAT expression", in which `nchar` is an integer variable. The error throwed out is like the following:

```txt
gfortran Error: Unexpected element ‘<’ in format string at (1)
```

To fix this issue, there are 2 common ways:

1. As GUN Fortran recommended, we can construct a format string dynamicly. See [here](http://gcc.gnu.org/onlinedocs/gfortran/Variable-FORMAT-expressions.html)

2. Simply replace `<nchar>` with a number which is large enough, say `10000`. Now, the `format` line will looks like this:

    ```fortran
    5     format(i2/i10/4e15.5,10000a1)
    ```

## Issue of character '$'

`gfortran` does not allow using `$` in variable name by default. Or it will throw the following error:

```txt
# when building xrayinvr
# in rayinvr/trc.f:498:17:
Fatal Error: Invalid character ‘$’ at (1). Use ‘-fdollar-ok’ to allow it as an extension compilation terminated
```

To fix this error, you should add option `-fdollar-ok` to `gfortran`. To be more detail, add this option to `FFLAGS` at the beginning of each `Makefile` file:

```make
FFLAGS = -O -fdollar-ok
```

## Change default foreground and backgroud color of Xwindow

`pltlib/xbuplot.c` uses `Xlib` to plot graphics. Foreground and background settings locates at line 99 and 100:

```c
  /*  create opaque window  */
  win = XCreateSimpleWindow (display, RootWindow (display, screen_num),
                x, y, win_width, win_height, border_width,
        WhitePixel (display, screen_num),    // foreground
        BlackPixel (display, screen_num));   // background
```

## Read v.in with higher precision

The default precision of `v.in` is `f7.2`. If you want to read a higher precision version of `v.in`(say, `f8.3`), just open the `main.f` that reads the `v.in`, change the lines:

```fortran
15       format(i2,1x,10f7.2)
25       format(3x,10i7)
25       format(3x,10(5x,i2))
```

to:

```fortran
15       format(i2,1x,10f8.3)
25       format(3x,10i8)
25       format(3x,10(6x,i2))
```

## Promoted X axis precision from 1m to 0.1m

Changed `rayinvr/main.f:890` (and many others):

```fortran
if(abs(xshotr-xshotf).lt..001.and.idr(is).eq.idf) then
```

to:

```fortran
if(abs(xshotr-xshotf).lt..0001.and.idr(is).eq.idf) then
```

## Some syntax errors

Changed `pltsyn/pltsec.f:253`:

```fortran
108     if(namp.eq.0) write(*,*, fmt="(/
```

to:

```fortran
108     if(namp.eq.0) write(*, fmt="(/
```

## Adaption for large models

Enlarged some parameters in `rayinvr/rayinvr.par` and `tramp/tramp.par` for processing larger models.

Promoted `prayi` to `100000` to enable huge `tx.in` files (>100000 lines).

## promote default `real` and `integer` precision from 4 bytes to 8 bytes

Add `-fdefault-real-8 -fdefault-double-8 -fdefault-integer-8` options to `FFLAGS` in `Makefile`, and replace fixed declaritions (`real*4`) in the code with flexible ones (`real`).
