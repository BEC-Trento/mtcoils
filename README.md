# mtcoils
Static B field calculator for arbitrary configurations of coils

## Yet another coil calculator
Exactly. To be honest, this started when I tried to re-use some piece of code written by a colleague of mine for the same task, and found it easier in the end to rewrite it from scratch instead than to compile it and make it work. At that point, most of the work was done and I thought is could be nice to finish this small piece of project. So here we are:


## Usage

    pip install --user git+https://github.com/BEC-Trento/mtcoils.git

Or simply clone the repo. Look at the [examples](examples) folder for some (commented) use case.

## The math

Ref. [math.pdf](doc/math.pdf)

## References and links

* [1] [Magnetic field of a coil in different coordinate systems](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20010038494.pdf)
* [Another example](http://docs.enthought.com/mayavi/mayavi/auto/example_magnetic_field.html#example-magnetic-field) for the same task, within the powerful [mayavi](http://docs.enthought.com/mayavi/mayavi/index.html) library
