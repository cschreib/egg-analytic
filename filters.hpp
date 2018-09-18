#ifndef FILTERS_INCLUDED
#define FILTERS_INCLUDED

#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

struct filter_options {
    std::string filter_db;
    bool filter_flambda = false;
    bool filter_photons = false;
    bool trim_filters = false;
};

class filter_database {
public :

    filter_database() = default;

    void read_options(program_arguments& opts) {
        filter_options fopts;
        opts.read(arg_list(
            fopts.filter_db, fopts.filter_flambda, fopts.filter_photons, fopts.trim_filters
        ));

        initialize(fopts);
    }

    void initialize(const filter_options& opts) {
        filter_db_file = opts.filter_db;
        filter_db = read_filter_db(filter_db_file);
        filter_flambda = opts.filter_flambda;
        filter_photons = opts.filter_photons;
        trim_filters = opts.trim_filters;
    }

    filter_t read_filter(const std::string& band) const {
        filter_t fil;
        vif_check(get_filter(filter_db, band, fil),
            "could not find filter '", band, "' in '", filter_db_file, "'");

        // Truncate
        if (trim_filters) {
            vec1u idg = where(fil.res/max(fil.res) > 1e-3);
            vif_check(!idg.empty(), "filter '", band, "' has no usable data");
            fil.lam = fil.lam[idg[0]-_-idg[-1]];
            fil.res = fil.res[idg[0]-_-idg[-1]];
        }

        // Apply filter definition
        if (filter_flambda) {
            // Filter is defined such that it must integrate f_lambda and not f_nu
            // f_lambda*r(lambda) ~ f_nu*[r(lambda)/lambda^2]
            fil.res /= sqr(fil.lam);
        }
        if (filter_photons) {
            // Filter is defined such that it integrates photons and not energy
            // n(lambda)*r(lambda) ~ f(lambda)*[r(lambda)*lambda]
            fil.res *= fil.lam;
        }

        // Re-normalize filter
        fil.res /= integrate(fil.lam, fil.res);
        fil.rlam = integrate(fil.lam, fil.lam*fil.res);

        return fil;
    }

private:

    bool filter_flambda = false;
    bool filter_photons = false;
    bool trim_filters = false;
    filter_db_t filter_db;
    std::string filter_db_file;
};

#endif
