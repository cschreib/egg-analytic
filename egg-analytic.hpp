#include <phypp.hpp>

namespace phypp {
namespace thread {
    template<typename W, typename T>
    struct worker_with_workspace {
        lock_free_queue<T> input;
        W                  wsp;
        std::atomic<bool>  shutdown;
        std::thread        impl;

        template<typename F, typename ... Args>
        explicit worker_with_workspace(const F& f, const Args&... args) : wsp(args...),
            shutdown(false), impl([this,f]() {

            T t;
            while (!shutdown) {
                while (input.pop(t)) {
                    f(wsp, t);
                }
            }
        }) {}

        ~worker_with_workspace() {
            join();
        }

        void join() {
            if (impl.joinable()) {
                impl.join();
            }
        }

        uint_t workload() const {
            return input.size();
        }
    };

    template<typename T>
    struct worker_no_workspace {
        lock_free_queue<T> input;
        std::atomic<bool>  shutdown;
        std::thread        impl;

        template<typename F>
        explicit worker_no_workspace(const F& f) : shutdown(false),
            impl([this,f]() {

            T t;
            while (!shutdown) {
                while (input.pop(t)) {
                    f(t);
                }
            }
        }) {}

        ~worker_no_workspace() {
            join();
        }

        void join() {
            if (impl.joinable()) {
                impl.join();
            }
        }

        uint_t workload() const {
            return input.size();
        }
    };

    template<typename T, typename W = void>
    struct worker_pool {
        using worker = typename std::conditional<std::is_same<W, void>::value,
            worker_no_workspace<T>, worker_with_workspace<W,T>>::type;

        std::vector<std::unique_ptr<worker>> workers;
        uint_t last_push = 0;

        worker_pool() = default;

        template<typename F, typename ... Args>
        explicit worker_pool(uint_t nthread, const F& f, const Args&... args) {
            start(nthread, f, args...);
        }

        ~worker_pool() {
            join();
        }

        template<typename F, typename ... Args>
        void start(uint_t nthread, const F& f, const Args&... args) {
            workers.clear();
            workers.reserve(nthread);
            for (uint_t i = 0; i < nthread; ++i) {
                workers.emplace_back(new worker(f, args...));
            }

            last_push = workers.size()-1;
        }

        void join() {
            for (uint_t i : range(workers)) {
                workers[i]->shutdown = true;
            }

            for (uint_t i : range(workers)) {
                workers[i]->join();
            }
        }

        void process(T t) {
            ++last_push;
            if (last_push >= workers.size()) {
                last_push = 0;
            }

            workers[last_push]->input.push(std::move(t));
        }

        void process(uint_t i, T t) {
            workers[i]->input.push(std::move(t));
        }

        void consume_all() {
            bool finished = false;
            while (!finished) {
                finished = true;
                for (uint_t i : range(workers)) {
                    if (workers[i]->workload() > 0) {
                        finished = false;
                        break;
                    }
                }

                sleep_for(1e-6);
            }
        }

        uint_t size() const {
            return workers.size();
        }

        uint_t remaining() const {
            uint_t nleft = 0;
            for (uint_t i : range(workers)) {
                nleft += workers[i]->workload();
            }

            return nleft;
        }
    };
}
}

namespace egg {
    namespace impl {
        struct csmf {
            vec1d zu, zl, zc;
            vec1d a_mstar1, a_mstar2, a_phistar1, a_phistar2, a_index1, a_index2;
            vec1d p_mstar1, p_mstar2, p_phistar1, p_phistar2, p_index1, p_index2;

            csmf() {
                // My stellar mass functions in GS
                zl = {0.3, 0.7, 1.2, 1.9, 2.5, 3.5};
                zu = {0.7, 1.2, 1.9, 2.5, 3.5, 4.5};

                // Active
                a_mstar1 = {11.0000, 11.0000, 11.0000, 11.0000, 11.0000, 11.0000};
                a_mstar2 = {10.6418, 10.7292, 10.6717, 10.8404, 10.9443, 11.0000};
                a_phistar1 = {0.000898887, 0.000718160, 0.000465684, 0.000213874, 0.000212404, 3.69e-05};
                a_index1 = {-1.40000, -1.40000, -1.50000, -1.57000, -1.60000, -1.80000};
                a_phistar2 = {8.30778e-05, 0.000404045, 0.000417749, 0.000406023, 9.06860e-05, 0.0};
                a_index2 = {0.500000, 0.500000, 0.500000, 0.00000, 0.500000, 0.500000};

                // Passive
                p_mstar1 = {11.0000, 11.0000, 11.0000, 11.0000, 11.0000, 11.0000};
                p_mstar2 = {11.0426, 10.8601, 10.8342, 11.0471, 10.9439, 11.0000};
                p_phistar1 = {7.77453e-05, 3.54586e-05, 2.29979e-05, 1.00000e-05, 0.00000, 0.00000};
                p_index1 = {-1.65000, -1.60000, -1.25000, -1.00000, -1.00000, -1.35000};
                p_phistar2 = {0.00154472, 0.00104263, 0.000624682, 0.000173119, 0.000122278, 3.00000e-05};
                p_index2 = {-0.481039, 0.0594024, 0.296244, -0.166611, -0.263124, -0.300000};

                // Baldry et al. z=0
                prepend(zu,         vec1d{min(zl)});
                prepend(zl,         vec1d{0.0});

                prepend(a_mstar1,   vec1d{10.72 + 0.25});
                prepend(a_mstar2,   vec1d{10.72 + 0.25});
                prepend(a_phistar1, vec1d{0.71e-3});
                prepend(a_phistar2, vec1d{0.0});
                prepend(a_index1,   vec1d{-1.45});
                prepend(a_index2,   vec1d{0.3});

                prepend(p_mstar1,   vec1d{10.72 + 0.25});
                prepend(p_mstar2,   vec1d{10.72 + 0.25});
                prepend(p_phistar1, vec1d{0.08e-3});
                prepend(p_phistar2, vec1d{3.25e-3});
                prepend(p_index1,   vec1d{-1.45});
                prepend(p_index2,   vec1d{-0.45});

                // Note: z > 4.0 is obtained by keeping the shape of the last redshift bin and
                //       decreasing phistar, following the total stellar mass density of
                //       Grazian et al. (2015)
                double zmax = 15.0;
                uint_t ilast = a_mstar1.size()-1;
                uint_t nhzb = ceil((zmax - max(zu))/0.2);
                double dhz = (zmax - max(zu))/nhzb;

                append(zl, dindgen(nhzb)*dhz + max(zu));
                append(zu, dindgen(nhzb)*dhz + max(zu) + dhz);

                auto g15_rhostar = vectorize_lambda([](double z){
                    return (z < 6 ? e10(-0.43*z) : e10(-0.43*6.0)*e10(-0.7*(z-6.0)));
                });

                vec1d decrease = g15_rhostar(0.5*(zl+zu)[(ilast+1)-_])/g15_rhostar(0.5*(zl[ilast] + zu[ilast]));
                vec1d index = interpolate({-1.8, -2.0, -2.1, -2.1}, {4.0, 5.0, 6.0, 6.01}, 0.5*(zl+zu)[(ilast+1)-_]);

                for (uint_t i : range(decrease)) {
                    a_mstar1.push_back(a_mstar1[ilast]);
                    a_mstar2.push_back(a_mstar2[ilast]);
                    a_index1.push_back(index[i]);
                    a_index2.push_back(a_index2[ilast]);
                    a_phistar1.push_back(decrease[i]*a_phistar1[ilast]);
                    a_phistar2.push_back(decrease[i]*a_phistar2[ilast]);

                    p_mstar1.push_back(p_mstar1[ilast]);
                    p_mstar2.push_back(p_mstar2[ilast]);
                    p_index1.push_back(p_index1[ilast]);
                    p_index2.push_back(p_index2[ilast]);
                    p_phistar1.push_back(decrease[i]*p_phistar1[ilast]);
                    p_phistar2.push_back(decrease[i]*p_phistar2[ilast]);
                }

                zc = 0.5*(zu + zl);
            }

            static double schechter2(double m, double mstar1, double index1, double phistar1,
                double mstar2, double index2, double phistar2) {
                double tm1 = e10(m-mstar1);
                double tm2 = e10(m-mstar2);
                return log(10.0)*(exp(-tm1)*phistar1*pow(tm1, 1+index1) + exp(-tm2)*phistar2*pow(tm2, 1+index2));
            }

            double evaluate_sf(double m, double z) {
                return schechter2(m,
                    interpolate(a_mstar1, zc, z),
                    interpolate(a_index1, zc, z),
                    max(interpolate(a_phistar1, zc, z), 0.0),
                    interpolate(a_mstar2, zc, z),
                    interpolate(a_index2, zc, z),
                    max(interpolate(a_phistar2, zc, z), 0.0)
                );
            }

            double evaluate_qu(double m, double z) {
                return schechter2(m,
                    interpolate(p_mstar1, zc, z),
                    interpolate(p_index1, zc, z),
                    max(interpolate(p_phistar1, zc, z), 0.0),
                    interpolate(p_mstar2, zc, z),
                    interpolate(p_index2, zc, z),
                    max(interpolate(p_phistar2, zc, z), 0.0)
                );
            }
        };
    }

    struct generator_options {
        // Main parameters of the survey
        std::string selection_band;
        double maglim = dnan;
        double area = 1.0; // [deg^2]
        vec1s filters;

        // Implementation parameters
        double logmass_min   = 4.0;
        double logmass_max   = 12.0;
        uint_t logmass_steps = 50;

        uint_t a_steps  = 50;
        uint_t b_steps  = 50;
        uint_t bt_steps = 5;

        uint_t seds_step = 1;

        bool naive_igm = false;

        // Resources
        std::string share_dir = "./";
        std::string filter_db = share_dir+"filter-db/db.dat";
        std::string sed_lib;
        std::string sed_lib_imf = "salpeter";
        bool filter_flambda = false;
        bool filter_photons = false;
        bool trim_filters = false;

        // Execution policy
        uint_t nthread = 0;
        uint_t max_queued_models = npos;
        bool strict_maglim = true;
    };

    class generator {
    public :
        bool initialized = false;

        // EGG SED library
        vec2b use;
        vec3f sed, lam;
        vec2f bvj, buv;
        vec3d flux;

        // Single SED library
        vec2f single_sed, single_lam;
        vec1b single_use;
        vec1u single_type;
        vec2d single_flux;
        vec1d single_mass;

        // Mass function
        impl::csmf cs_mf;

        // Internal variables
        double flim = dnan;
        cosmo_t cosmo;
        double area = dnan;

        // Config
        filter_db_t filter_db;
        std::string selection_band;
        filter_t selection_filter;
        vec1s bands;
        vec<1,filter_t> filters;
        bool naive_igm = false;
        bool filter_flambda = false;
        bool filter_photons = false;
        bool trim_filters = false;
        uint_t nthread = 0;
        uint_t max_queued_models = npos;
        bool strict_maglim = true;
        bool single_sed_library = false;

        // Base arrays
        vec1d m, nm;
        double dm;
        vec1d a, b;
        vec1d uv, vj;
        vec1d bt;
        uint_t seds_step = 1;

        // Base distributions
        static double gaussian_integrate(double x1, double x2, double mu, double sigma) {
            // Sample (incorect!)
            // return exp(-sqr(x-mu)/(2.0*sqr(sigma)))/(sigma*sqrt(2*dpi));

            // Integrate
            return 0.5*(erf((x2 - mu)/(sigma*sqrt(2.0))) - erf((x1 - mu)/(sigma*sqrt(2.0))));
        }

        static vec1d gaussian(const vec1d& x, double mu, double sigma) {
            vec1d ret(x.dims);
            for (uint_t ix : range(x)) {
                const double x1 = (ix == 0          ? x.safe[ix] : 0.5*(x.safe[ix] + x.safe[ix-1]));
                const double x2 = (ix == x.size()-1 ? x.safe[ix] : 0.5*(x.safe[ix] + x.safe[ix+1]));
                ret.safe[ix] = gaussian_integrate(x1, x2, mu, sigma);
            }

            return ret;
        }

        static double lognormal_integrate(double x1, double x2, double mu, double sigma) {
            // Sample (incorect!)
            // return exp(-sqr(log(x/mu))/(2.0*sqr(sigma)))/(x*sigma*sqrt(2*dpi));

            // Integrate
            double up = erf(log(x2/mu)/(sigma*sqrt(2.0)));
            double low = (x1 > 0.0 ? erf(log(x1/mu)/(sqrt(2.0)*sigma)) : -1.0);
            return 0.5*(up - low);
        }

        static vec1d lognormal(const vec1d& x, double mu, double sigma) {
            vec1d ret(x.dims);
            for (uint_t ix : range(x)) {
                const double x1 = (ix == 0          ? x.safe[ix] : 0.5*(x.safe[ix] + x.safe[ix-1]));
                const double x2 = (ix == x.size()-1 ? x.safe[ix] : 0.5*(x.safe[ix] + x.safe[ix+1]));
                ret.safe[ix] = lognormal_integrate(x1, x2, mu, sigma);
            }

            return ret;
        }

        static double upper_lognormal(double xlow, double mu, double sigma) {
            return 0.5 - 0.5*erf(log(xlow/mu)/(sigma*sqrt(2.0)));
        }

        static double upper_gaussian(double xlow, double mu, double sigma) {
            return 0.5 - 0.5*erf((xlow - mu)/(sigma*sqrt(2.0)));
        }

        static double lower_gaussian(double xup, double mu, double sigma) {
            return 0.5 + 0.5*erf((xup - mu)/(sigma*sqrt(2.0)));
        }

        generator() = default;
        virtual ~generator() = default;

        bool read_filter(const std::string& band, filter_t& fil) const {
            if (!get_filter(filter_db, band, fil)) {
                return false;
            }

            // Truncate
            if (trim_filters) {
                vec1u idg = where(fil.res/max(fil.res) > 1e-3);
                phypp_check(!idg.empty(), "filter '", band, "' has no usable data");
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

            return true;
        }

        void initialize(const generator_options& opts) {
            phypp_check(!opts.selection_band.empty(),
                "please provide the name of the selection band in the options");
            phypp_check(is_finite(opts.maglim),
                "please provide the value for the magnitude limit in the options");

            // Read SED library
            naive_igm = opts.naive_igm;
            std::string sed_file;
            if (opts.sed_lib.empty()) {
                sed_file = opts.share_dir+"opt_lib_fast_hd_noigm.fits";
                if (naive_igm) {
                    sed_file = opts.share_dir+"opt_lib_fast_hd.fits";
                }
            } else {
                sed_file = opts.sed_lib;
            }

            {
                fits::input_table itbl(sed_file);

                fits::column_info cinfo;
                phypp_check(itbl.read_column_info("sed", cinfo),
                    "invalid SED library, must have SED column");
                phypp_check(cinfo.dims.size() == 2 || cinfo.dims.size() == 3,
                    "invalid SED library, SED column must have 2 or 3 dimensions");

                if (cinfo.dims.size() == 2) {
                    note("found a single SED library (no UVJ color binning)");
                    single_sed_library = true;
                    itbl.read_columns("use", single_use, "lam", single_lam,
                        "sed", single_sed, "type", single_type);

                    if (itbl.read_column("mass", single_mass)) {
                        note("SED library has a mass calibration");
                    }

                    buv.resize(2,0);
                    bvj.resize(2,0);
                } else {
                    note("found an EGG library (with UVJ color binning)");
                    single_sed_library = false;
                    itbl.read_columns(ftable(use, lam, sed, buv, bvj));
                }
            }

            // The SEDs in the library are given in Lsun/Msun.
            // Therefore they assume an IMF. We need to correct for this
            // if this assumed IMF is different from the one used to
            // build the stellar mass functions (i.e., Salpeter).
            // Here we assume this is a simple global rescaling factor.
            double conv = 1.0;
            if (opts.sed_lib_imf == "salpeter") {
                conv = 1.0; // OK, EGG assumes Salpeter throughout
            } else if (opts.sed_lib_imf == "chabrier") {
                conv = e10(-0.2);
            }

            if (single_sed_library) {
                single_sed *= conv;
            } else {
                sed *= conv;
            }

            // Skip SEDs (if asked)
            seds_step = opts.seds_step;
            if (single_sed_library) {
                for (uint_t ised : range(single_use)) {
                    if (ised % seds_step != 0) {
                        single_use.safe[ised] = false;
                    }
                }

                note("generating mock with ", count(single_use), " SEDs");
            } else {
                for (uint_t iuv : range(use.dims[0]))
                for (uint_t ivj : range(use.dims[1])) {
                    if ((iuv+ivj) % seds_step != 0) {
                        use.safe(iuv, ivj) = false;
                    }
                }

                note("generating mock with ", count(use), " disk and bulge SEDs");
            }

            // Internal parameters
            flim = e10(0.4*(23.9 - opts.maglim));
            cosmo = get_cosmo("std"); // must be specific cosmology for which EGG is calibrated
            area = opts.area*sqr(dpi/180.0)/(4*dpi); // fraction of the sky we simulate

            // Base arrays
            m  = rgen(opts.logmass_min, opts.logmass_max, opts.logmass_steps);
            dm = m[1] - m[0];
            nm = e10(m);

            a  = rgen(-5.0, 6.0, opts.a_steps);
            b  = rgen(-0.1, 0.2, opts.b_steps);

            uv = bin_center(buv);
            vj = bin_center(bvj);

            bt = rgen(0.0, 1.0, opts.bt_steps);

            // Read filter library
            filter_db = read_filter_db(opts.filter_db);
            filter_flambda = opts.filter_flambda;
            filter_photons = opts.filter_photons;
            trim_filters = opts.trim_filters;

            // Find selection filter
            selection_band = opts.selection_band;
            phypp_check(read_filter(selection_band, selection_filter),
                "could not find selection filter, aborting");

            // Find flux filters
            bands = opts.filters;
            filters.resize(bands.size());
            for (uint_t l : range(bands)) {
                phypp_check(read_filter(bands[l], filters[l]),
                    "could not find filter '", bands[l], "', aborting");
            }

            // Other parameters
            nthread = opts.nthread;
            max_queued_models = opts.max_queued_models;
            if (max_queued_models == npos) {
                // Automatic: keep at most two queued models per thread
                max_queued_models = 2*nthread;
            }

            strict_maglim = opts.strict_maglim;

            initialized = true;
        }

        virtual void on_disk_sed_changed(uint_t ised_d) {}
        virtual void on_bulge_sed_changed(uint_t ised_b) {}
        virtual void on_generated(uint_t iter, uint_t im, uint_t it, uint_t ised_d, uint_t ised_b,
            uint_t ibt, double ngal, const vec1d& fdisk, const vec1d& fbulge) = 0;

        void apply_madau_igm(double z, const vec1d& wl, vec1d& flx) {
            // http://adsabs.harvard.edu/abs/1995ApJ...441...18M
            // TODO: check this implementation someday, I suspect this is wrong or
            // very approximate (taken directly from FAST)

            double da; {
                double l0 = 1050.0*(1.0 + z);
                double l1 = 1170.0*(1.0 + z);
                uint_t nstep = 100;
                vec1d tl = rgen(l0, l1, nstep);
                vec1d ptau = exp(-3.6e-3*pow(tl/1216.0, 3.46));
                da = mean(ptau);
            }

            double db; {
                double l0 = 920.0*(1.0 + z);
                double l1 = 1015.0*(1.0 + z);
                uint_t nstep = 100;
                vec1d tl = rgen(l0, l1, nstep);
                vec1d ptau = exp(-1.7e-3*pow(tl/1026.0, 3.46) - 1.2e-3*pow(tl/972.5, 3.46) -
                    9.3e-4*pow(tl/950.0, 3.46));
                db = mean(ptau);
            }

            uint_t l0 = lower_bound(wl, 0.0912);
            uint_t l1 = lower_bound(wl, 0.1026);
            uint_t l2 = lower_bound(wl, 0.1216);

            for (uint_t l : range(l0)) {
                flx.safe[l] = 0;
            }
            for (uint_t l : range(l0, l1)) {
                flx.safe[l] *= db;
            }
            for (uint_t l : range(l1, l2)) {
                flx.safe[l] *= da;
            }
        }

        struct generated_data {
            uint_t iter, im, it, ised_d, ised_b, ibt;
            double ngal;
            vec1d fdisk, fbulge;
        };

        void generate(double zf, double dz) {
            phypp_check(initialized, "please first call initialize() with your desired survey setup");

            // Pre-compute stuff
            const vec1d abar = [&]() {
                const vec1d a0 = 0.58*erf(m - 10.0) + 1.39;
                const vec1d a1 = -0.34 + 0.3*max(m - 10.35, 0.0);
                return min(a0 + a1*min(zf, 3.3), 2.0);
            }();
            const vec1d asig = 0.1 + 0.3*clamp(zf-1.0, 0, 1)*(0.17 + clamp(1.0 - abs(m - 10.3), 0, 1));
            const vec1d bbar = 0.1*(m - 11.0);

            const double col_cor = 0.2*(zf < 0.5 ? (0.5 - zf)/0.5 : 0.0);
            const double ML_cor = e10(-interpolate({0.15,0.15,0.0,0.0,-0.6}, {0.0,0.45,1.3,6.0,8.0}, zf));

            const double dvdz = (vuniverse(zf+dz/2.0, cosmo) - vuniverse(zf-dz/2.0, cosmo))*(area/dz);

            double df = lumdist(zf, cosmo);

            // Pre-compute fluxes in selection band in library
            uint_t first_sed = npos;
            double fmax = 0.0;
            if (single_sed_library) {
                single_flux.resize(single_use.dims, filters.size()+1);
                for (uint_t ised : range(single_use)) {
                    if (!single_use.safe[ised]) continue;

                    if (first_sed == npos) {
                        first_sed = ised;
                    }

                    vec1d tlam = single_lam.safe(ised,_);
                    vec1d tsed = single_sed.safe(ised,_);

                    if (!naive_igm) {
                        apply_madau_igm(zf, tlam, tsed);
                    }

                    tsed = lsun2uJy(zf, df, tlam, tsed);
                    tlam *= (1.0 + zf);

                    single_flux.safe(ised,0) = ML_cor*sed2flux(
                        selection_filter.lam, selection_filter.res, tlam, tsed
                    );

                    if (single_flux.safe(ised,0) > fmax) {
                        fmax = single_flux.safe(ised,0);
                    }

                    for (uint_t l : range(filters)) {
                        single_flux.safe(ised,l+1) = ML_cor*sed2flux(
                            filters[l].lam, filters[l].res, tlam, tsed
                        );
                    }
                }
            } else {
                flux.resize(use.dims, filters.size()+1);
                for (uint_t iuv : range(use.dims[0]))
                for (uint_t ivj : range(use.dims[1])) {
                    if (!use.safe(iuv,ivj)) continue;

                    if (first_sed == npos) {
                        first_sed = iuv*use.dims[1] + ivj;
                    }

                    vec1d tlam = lam.safe(iuv,ivj,_);
                    vec1d tsed = sed.safe(iuv,ivj,_);

                    if (!naive_igm) {
                        apply_madau_igm(zf, tlam, tsed);
                    }

                    tsed = lsun2uJy(zf, df, tlam, tsed);
                    tlam *= (1.0 + zf);

                    flux.safe(iuv,ivj,0) = ML_cor*sed2flux(
                        selection_filter.lam, selection_filter.res, tlam, tsed
                    );

                    if (flux.safe(iuv,ivj,0) > fmax) {
                        fmax = flux.safe(iuv,ivj,0);
                    }

                    for (uint_t l : range(filters)) {
                        flux.safe(iuv,ivj,l+1) = ML_cor*sed2flux(
                            filters[l].lam, filters[l].res, tlam, tsed
                        );
                    }
                }
            }

            // Select probability threshold below which SEDs are ignored
            const uint_t nsed = (single_sed_library ? count(single_use) : count(use));
            const double pthr = 0.1/nsed;

            // Find minimum mass we need to bother with
            const double mmin = log10(flim/fmax);
            const uint_t m0 = lower_bound(m, mmin);

            thread::worker_pool<generated_data> pool;
            if (nthread > 0) {
                // Start fitting threads
                pool.start(nthread, [this](const generated_data& d) {
                    on_generated(
                        d.iter, d.im, d.it, d.ised_d, d.ised_b, d.ibt, d.ngal, d.fdisk, d.fbulge
                    );
                });
            }

            uint_t iter = 0;

            // Processing function (dispatch to direct processing or thread pool)
            auto process = [&](uint_t tim, uint_t tit, uint_t tised_d, uint_t tised_b,
                uint_t tibt, double tngal, const vec1d& tfdisk, const vec1d& tfbulge) {

                if (strict_maglim && tfdisk.safe[0] + tfbulge.safe[0] < flim) {
                    // Do not send galaxies below the magnitude limit
                    return;
                }

                if (nthread == 0) {
                    // No multi-threading, handle now
                    on_generated(iter, tim, tit, tised_d, tised_b, tibt, tngal, tfdisk, tfbulge);
                } else {
                    // Multi-threading.
                    // First make sure we don't have too many queued calculations,
                    // and if so wait for a moment until previous calculations have finished.
                    while (max_queued_models > 0 && pool.remaining() > max_queued_models) {
                        thread::sleep_for(1e-6);
                    }

                    // Now send the model to the thread pool
                    pool.process(generated_data{iter, tim, tit, tised_d, tised_b,
                        tibt, tngal, tfdisk, tfbulge});
                }

                ++iter;
            };

            // Integrate over stellar masses
            for (uint_t im : range(m0, m.size())) {
                const double mm = nm.safe[im];

                if (single_sed_library) {
                    // Simpler approach where we only generate fluxes from the SED library
                    // as it is, not trying to add bulge and disk components

                    // Integrate over galaxy type (SF=1 or QU=0)
                    for (uint_t it : {0, 1}) {
                        // Probability
                        const double nmz = dz*dm*dvdz*(it == 0 ?
                            cs_mf.evaluate_qu(m.safe[im], zf) : cs_mf.evaluate_sf(m.safe[im], zf));

                        if (!single_mass.empty()) {
                            // We have a mass calibration, pick the SEDs with masses within our bin
                            vec1u imass = where(abs(single_mass - m.safe[im]) < 0.5*dm);
                            if (imass.empty()) {
                                // No SED with this mass, pick the closest one
                                uint_t iclose = min_id(abs(m.safe[im] - single_mass));
                                imass = where(abs(single_mass - single_mass[iclose]) < 0.5*dm);
                            }

                            for (uint_t ised : imass) {
                                // Notify disk SED has changed
                                if (nthread == 0) {
                                    on_disk_sed_changed(ised);
                                }

                                // Load flux
                                const vec1d fdisk = mm*single_flux.safe(ised,_);

                                // Forward to derived class
                                process(im, it, ised, first_sed, 0, nmz, fdisk, fdisk*0.0);
                            }
                        } else {
                            for (uint_t ised : range(single_use)) {
                                if (!single_use.safe[ised] || single_type.safe[ised] != it) continue;

                                // Notify disk SED has changed
                                if (nthread == 0) {
                                    on_disk_sed_changed(ised);
                                }

                                // Load flux
                                const vec1d fdisk = mm*single_flux.safe(ised,_);

                                // Forward to derived class
                                process(im, it, ised, first_sed, 0, nmz/nsed, fdisk, fdisk*0.0);
                            }
                        }
                    }
                } else {
                    // EGG approach

                    // pblue
                    vec2d pblue(use.dims); {
                        vec1d pa = gaussian(a, abar.safe[im], asig.safe[im]);
                        for (uint_t ia : range(a)) {
                            vec1d puv = gaussian(uv, 2*col_cor + 0.45 + 0.545*a.safe[ia], 0.12);
                            vec1d pvj = gaussian(vj,   col_cor +        0.838*a.safe[ia], 0.12);
                            for (uint_t iuv : range(uv))
                            for (uint_t ivj : range(vj)) {
                                pblue.safe(iuv,ivj) += puv.safe[iuv]*pvj.safe[ivj]*pa.safe[ia];
                            }
                        }

                        pblue[where(!use)] = 0;
                        pblue /= total(pblue);
                    }

                    // pred
                    vec2d pred(use.dims); {
                        vec1d pb = gaussian(b, bbar.safe[im], 0.1);
                        pb.safe[0]           += lower_gaussian(-0.1, bbar.safe[im], 0.10);
                        pb.safe[pb.size()-1] += upper_gaussian( 0.2, bbar.safe[im], 0.10);
                        for (uint_t ib : range(b)) {
                            vec1d puv = gaussian(uv, 2*col_cor + 1.85 + 0.88*b.safe[ib], 0.10);
                            vec1d pvj = gaussian(vj,   col_cor + 1.25 +      b.safe[ib], 0.10);
                            for (uint_t iuv : range(uv))
                            for (uint_t ivj : range(vj)) {
                                pred.safe(iuv,ivj) += puv.safe[iuv]*pvj.safe[ivj]*pb.safe[ib];
                            }
                        }

                        pred[where(!use)] = 0;
                        pred /= total(pred);
                    }

                    // Integrate over galaxy type (SF or QU)
                    for (uint_t it : {0, 1}) {
                        // N(M*,t,z)
                        const double nmz = dz*dm*dvdz*(it == 0 ?
                            cs_mf.evaluate_qu(m.safe[im], zf) : cs_mf.evaluate_sf(m.safe[im], zf));

                        // p(B/T | M*,t)
                        const double btbar = (it == 0 ? 0.5*pow(mm/1e10, 0.10) : 0.2*pow(mm/1e10, 0.27));
                        vec1d pbt = lognormal(bt, btbar, 0.2);
                        pbt.safe[bt.size()-1] += upper_lognormal(1.0, btbar, 0.2);

                        // Integrate over disk SED
                        for (uint_t ised_d : range(use)) {
                            if (!use.safe[ised_d]) continue;

                            const uint_t iduv = ised_d / use.dims[1];
                            const uint_t idvj = ised_d % use.dims[1];

                            // Notify disk SED has changed
                            if (nthread == 0) {
                                on_disk_sed_changed(ised_d);
                            }

                            // pdisk(d | z,M*) * N(M*,t,z)
                            const double pdisk = nmz*pblue.safe(iduv, idvj);
                            const double pred_bt1 = nmz*pred.safe(iduv, idvj);

                            // Load stuff
                            const vec1d fdisk = mm*flux.safe(iduv,idvj,_);

                            // Optimization: send B/T=0 now, don't need to iterate over SED bulge
                            if (pdisk >= nmz*pthr) {
                                process(im, it, ised_d, first_sed,
                                    0, pbt.safe[0]*pdisk, fdisk, fdisk*0.0);
                            }

                            // Optimization: send B/T=1 now, don't need to iterate over SED disk
                            // (trick: pretend we're actually iterating over SED bulge now)
                            if (pred_bt1 >= nmz*pthr) {
                                if (nthread == 0) {
                                    on_bulge_sed_changed(ised_d);
                                }

                                process(im, it, first_sed, ised_d,
                                    bt.size()-1, pbt.safe[bt.size()-1]*pred_bt1, fdisk*0.0, fdisk);
                            }

                            // Skip improbable SEDs
                            if (pdisk < nmz*pthr) continue;

                            // Integrate over bulge SED
                            for (uint_t ised_b : range(use)) {
                                if (!use.safe[ised_b]) continue;

                                const uint_t ibuv = ised_b / use.dims[1];
                                const uint_t ibvj = ised_b % use.dims[1];

                                // Notify bulge SED has changed
                                if (nthread == 0) {
                                    on_bulge_sed_changed(ised_b);
                                }

                                // p(blue) * pdisk(d | z,M*) * N(M*,t,z)
                                const double pb = pdisk*pblue.safe(ibuv,ibvj);
                                const double pr = pdisk*pred.safe(ibuv,ibvj);

                                // Skip improbable SEDs
                                if (pb < pdisk*pthr && pr < pdisk*pthr) continue;

                                // Load stuff
                                const vec1d fbulge = mm*flux.safe(ibuv,ibvj,_);

                                // Integrate over B/T (skipping B/T=0 and 1, already done above)
                                for (uint_t ibt : range(1, bt.size()-1)) {
                                    const double btn = bt.safe[ibt];
                                    const double bti = 1.0 - btn;

                                    // pbulge(b | z,M*,B/T)
                                    const double pbulge = (btn >= 0.6 ? pr : 0.5*(pr + pb));

                                    // Compute this galaxy's number density (per delta z)
                                    const double tprob = pbulge*pbt.safe[ibt];

                                    // And forward to derived class.
                                    process(im, it, ised_d, ised_b, ibt, tprob, fdisk*bti, fbulge*btn);
                                }
                            }
                        }
                    }
                }
            }
        }
    };
}
