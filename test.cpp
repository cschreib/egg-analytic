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

    // Generate population at the specified redshift
    gen.generate(z, 0.01);

    // We have our result!
    print("the average log stellar mass is: ", gen.mass_tot/gen.ngal);

    return 0;
}
