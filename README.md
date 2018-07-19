# EGG (analytic version)

This is an analytic version of the Empirical Galaxy Generator (EGG):

https://github.com/cschreib/egg

It is provided as a header-only C++ library, and requires some resources from the original EGG.


# Installing

Because this is header-only library, it does not need to be installed. However it has some dependencies which you do need to install before using it. In particular you must have the phy++ library installed:

https://github.com/cschreib/phypp

You must also download the original version of EGG (installing it is not required, we just need some files from there):

https://github.com/cschreib/egg


# Simple example of usage

Create a new file called ``test.cpp`` and copy the following:

```c++
#include "egg-analytic.hpp"

// Our generator class
class my_generator : public egg::generator {
public :
    // Accumulating things
    double ngal = 0, mass_tot = 0;

    my_generator() : egg::generator() {}

    void on_generated(uint_t im, uint_t it, uint_t ised_d, uint_t ised_b,
        uint_t ibt, double tngal, double fdisk, double fbulge) override {

        if (ftot >= flim) {
            // We passed the magnitude cut!

            // Sum things up
            ngal += tngal;
            mass_tot += m.safe[im]*tngal;
        }
    }
};

int phypp_main(int argc, char* argv[]) {
    // Configurable options
    double maglim = 24.5;
    std::string selection_band = "euclid-vis";
    double z = 1.0;
    std::string egg_dir = "./";

    // Read them from command line
    read_args(argc, argv, arg_list(selection_band, maglim, z, egg_dir));

    // Create our generator
    my_generator gen;

    // Setup the survey
    egg::generator_options opts;
    opts.share_dir = file::directorize(egg_dir);
    opts.filter_db = opts.share_dir+"filter-db/db.dat";
    opts.selection_band = selection_band;
    opts.maglim = maglim;
    gen.initialize(opts);

    // Generate population at the specified redshift
    gen.generate(z, 0.01);

    // We have our result!
    print("the average log stellar mass is: ", gen.mass_tot/gen.ngal);

    return 0;
}
```

Copy the ``egg-analytic.hpp`` header in the same directory, then run the following:

```bash
cphy++ optimize test.cpp
```

Now the program is compiled. To run it, do:

```bash
./test selection_band=hst-f160w maglim=26 z=2.5 egg_dir=<path/to/egg/share>
```

After a few seconds, this should print:

```bash
the average log stellar mass is: 9.53505
```

In the code above, the most important part is the ``on_generated()`` function. It is called at each step of the integration, for each value of the stellar mass, galaxy type (SF/QU), disk SED, bulge SED, and bulge-to-total ratio. It provides the number density of objects with these properties (``tngal``, which is mathematically ``dN/dz``), the flux of the disk component (``fdisk``) and the flux of the bulge component (``fbulge``).

You can access the corresponding physical quantites as follows:

* log stellar mass: ``m.safe[im]``.
* galaxy type: ``it`` (0: quiescent, 1: star-forming).
* disk SED: ``ised_d`` (index in the SED library).
* bulge SED: ``ised_b`` (index in the SED library).
* bulge-to-total ratio: ``bt.safe[ibt]``.

Note that you have to check yourself that the predicted flux for these galaxies is higher than the chosen magnitude cut.


# Going further

The example above is very simplistic, and simply computes the average log stellar mass at a given redshift, given a magnitude cut in a given band. This is not particularly useful. Here we will improve this example and instead compute the lowest stellar mass for which the sample is 90% complete.

For this, we will need to store, for each stellar mass in the grid, two values: the total number of *expected* galaxies, and the total number of *detected* galaxies. We will therefore create two arrays:

```c++
#include "egg-analytic.hpp"

// Our generator class
class my_generator : public egg::generator {
public :
    // Accumulating things
    vec1d ntot, ndet;        // <<---- ADDED

    // [...]
};
```

We then initialize them in the main function:

```c++
int phypp_main(int argc, char* argv[]) {
    // [...]

    // Create our generator
    my_generator gen;

    // Setup the survey
    egg::generator_options opts;
    opts.share_dir = file::directorize(egg_dir);
    opts.filter_db = opts.share_dir+"filter-db/db.dat";
    opts.selection_band = selection_band;
    opts.maglim = maglim;
    gen.initialize(opts);

    // Initialize our new arrays
    // Must be done after initialize(), because this is where the
    // stellar mass array is created.
    gen.ntot.resize(gen.m.size());  // <<---- ADDED
    gen.ndet.resize(gen.m.size());  // <<---- ADDED

    // [...]
}
```

And we add the code to locate the lowest mass with 80% completeness:

```c++
int phypp_main(int argc, char* argv[]) {
    // [...]

    // Generate population at the specified redshift
    gen.generate(z, 0.01);

    // Now we need to find the lowest mass at which the completeness is >90%
    uint_t i80 = gen.m.size();
    while (i80 > 0) {
        --i80;

        if (gen.ndet[i80]/gen.ntot[i80] < 0.9) {
            ++i80;
            break;
        }
    }

    // We have our result!
    print("the stellar mass 80% completeness is: ", gen.m[i80]);

    return 0;
}
```

We're left with the most important thing: the ``on_generated()`` function:

```c++
    void on_generated(uint_t im, uint_t it, uint_t ised_d, uint_t ised_b,
        uint_t ibt, double tngal, double fdisk, double fbulge) override {

        // Sum up all galaxies in 'ntot'
        ntot.safe[im] += tngal;

        // Compute total flux
        double ftot = fdisk + fbulge;
        if (ftot >= flim) {
            // We made it above the detection limit!

            // Sum up the detected galaxies in 'ndet'
            ndet.safe[im] += tngal;
        }
    }
```

Here is the final file:

```c++
#include "egg-analytic.hpp"

// Our generator class
class my_generator : public egg::generator {
public :
    // Accumulating things
    vec1d ntot, ndet;

    my_generator() : egg::generator() {}

    void on_generated(uint_t im, uint_t it, uint_t ised_d, uint_t ised_b,
        uint_t ibt, double tngal, double fdisk, double fbulge) override {

        // Sum up all galaxies in 'ntot'
        ntot.safe[im] += tngal;

        double ftot = fdisk + fbulge;
        if (ftot >= flim) {
            // Sum up only the detected galaxies in 'ndet'
            ndet.safe[im] += tngal;
        }
    }
};

int phypp_main(int argc, char* argv[]) {
    // Configurable options
    double maglim = 24.5;
    std::string selection_band = "euclid-vis";
    double z = 1.0;
    std::string egg_dir = "./";

    // Read them from command line
    read_args(argc, argv, arg_list(selection_band, maglim, z, egg_dir));

    // Create our generator
    my_generator gen;

    // Setup the survey
    egg::generator_options opts;
    opts.share_dir = file::directorize(egg_dir);
    opts.filter_db = opts.share_dir+"filter-db/db.dat";
    opts.selection_band = selection_band;
    opts.maglim = maglim;
    gen.initialize(opts);

    // Initialize our new arrays
    // Must be done after initialize(), because this is where the
    // stellar mass array is created.
    gen.ntot.resize(gen.m.size());
    gen.ndet.resize(gen.m.size());

    // Generate population at the specified redshift
    gen.generate(z, 0.01);

    // Now we need to find the lowest mass at which the completeness is >80%
    uint_t i80 = gen.m.size();
    while (i80 > 0) {
        --i80;

        if (gen.ndet[i80]/gen.ntot[i80] < 0.9) {
            ++i80;
            break;
        }
    }

    // We have our result!
    print("the stellar mass 90% completeness is: ", gen.m[i80]);

    return 0;
}
```

Compiling and running gives:

```bash
bash> ./test2 selection_band=hst-f160w maglim=26 z=2.2 egg_dir=<path/to/egg/share>
the stellar mass 90% completeness is: 9.18182
```
