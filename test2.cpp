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
    std::string selection_band = "subaru-B";
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
    uint_t i90 = gen.m.size();
    while (i90 > 0) {
        --i90;

        if (gen.ndet[i90]/gen.ntot[i90] < 0.9) {
            ++i90;
            break;
        }
    }

    // We have our result!
    print("the stellar mass 90% completeness is: ", gen.m[i90]);

    return 0;
}
