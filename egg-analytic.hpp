#include <phypp.hpp>

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
                prepend(p_phistar1, vec1d{3.25e-3});
                prepend(p_phistar2, vec1d{0.08e-3});
                prepend(p_index1,   vec1d{-0.45});
                prepend(p_index2,   vec1d{-1.45});

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
        double logmass_max   = 13.0;
        uint_t logmass_steps = 50;

        uint_t a_steps  = 50;
        uint_t b_steps  = 50;
        uint_t bt_steps = 5;

        uint_t seds_step = 1;

        // Resources
        std::string share_dir = "./";
        std::string filter_db = share_dir+"filter-db/db.dat";
        bool filter_flambda = false;
        bool filter_photons = false;
        bool trim_filters = false;
    };

    class generator {
    public :
        bool initialized = false;

        // SED library
        vec2b use;
        vec3f sed, lam;
        vec2f bvj, buv;
        vec3d flux;

        // Mass function
        impl::csmf cs_mf;

        // Internal variables
        double flim = dnan;
        cosmo_t cosmo;
        double area = dnan;
        std::string selection_band;
        filter_t selection_filter;
        vec1s bands;
        vec<1,filter_t> filters;

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

        void initialize(const generator_options& opts) {
            phypp_check(!opts.selection_band.empty(),
                "please provide the name of the selection band in the options");
            phypp_check(is_finite(opts.maglim),
                "please provide the value for the magnitude limit in the options");

            // Read SED library
            fits::read_table(opts.share_dir+"opt_lib_fast_hd.fits",
                ftable(use, lam, sed, buv, bvj));

            // Internal parameters
            flim = e10(0.4*(23.9 - opts.maglim));
            cosmo = get_cosmo("std"); // must be specific cosmology for which EGG is calibrated
            area = opts.area*sqr(dpi/180.0)/(4*dpi); // fraction of the sky we simulate

            // Base arrays
            seds_step = opts.seds_step;
            m  = rgen(opts.logmass_min, opts.logmass_max, opts.logmass_steps);
            dm = m[1] - m[0];
            nm = e10(m);

            a  = rgen(-5.0, 6.0, opts.a_steps);
            b  = rgen(-0.1, 0.2, opts.b_steps);

            uv = bin_center(buv);
            vj = bin_center(bvj);

            bt = rgen(0.0, 1.0, opts.bt_steps);

            // Read filter library
            auto filter_db = read_filter_db(opts.filter_db);

            // Find selection filter
            selection_band = opts.selection_band;
            phypp_check(get_filter(filter_db, selection_band, selection_filter),
                "could not find selection filter, aborting");

            bands = opts.filters;
            phypp_check(get_filters(filter_db, bands, filters),
                "could not find some filters, aborting");

            // Truncate filters
            if (opts.trim_filters) {
                vec1u idg = where(selection_filter.res/max(selection_filter.res) > 1e-3);
                phypp_check(!idg.empty(), "filter '", opts.selection_band, "' has no usable data");
                selection_filter.lam = selection_filter.lam[idg[0]-_-idg[-1]];
                selection_filter.res = selection_filter.res[idg[0]-_-idg[-1]];
                for (uint_t l : range(filters)) {
                    idg = where(filters[l].res/max(filters[l].res) > 1e-3);
                    phypp_check(!idg.empty(), "filter '", opts.filters[l], "' has no usable data");
                    filters[l].lam = filters[l].lam[idg[0]-_-idg[-1]];
                    filters[l].res = filters[l].res[idg[0]-_-idg[-1]];
                }
            }

            // Apply filter definition
            if (opts.filter_flambda) {
                // Filter is defined such that it must integrate f_lambda and not f_nu
                // f_lambda*r(lambda) ~ f_nu*[r(lambda)/lambda^2]
                selection_filter.res /= sqr(selection_filter.lam);
                for (uint_t l : range(filters)) {
                    filters[l].res /= sqr(filters[l].lam);
                }
            }
            if (opts.filter_photons) {
                // Filter is defined such that it integrates photons and not energy
                // n(lambda)*r(lambda) ~ f(lambda)*[r(lambda)*lambda]
                selection_filter.res *= selection_filter.lam;
                for (uint_t l : range(filters)) {
                    filters[l].res *= filters[l].lam;
                }
            }

            // Re-normalize filter
            selection_filter.res /= integrate(selection_filter.lam, selection_filter.res);
            selection_filter.rlam = integrate(selection_filter.lam, selection_filter.lam*selection_filter.res);
            for (uint_t l : range(filters)) {
                filters[l].res /= integrate(filters[l].lam, filters[l].res);
                filters[l].rlam = integrate(filters[l].lam, filters[l].res*filters[l].lam);
            }

            initialized = true;
        }

        virtual void on_disk_sed_changed(uint_t ised_d) {}
        virtual void on_bulge_sed_changed(uint_t ised_b) {}
        virtual void on_generated(uint_t im, uint_t it, uint_t ised_d, uint_t ised_b,
            uint_t ibt, double ngal, const vec1d& fdisk, const vec1d& fbulge) = 0;

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
            flux.resize(use.dims, filters.size()+1);
            for (uint_t iuv : range(flux.dims[0]))
            for (uint_t ivj : range(flux.dims[1])) {
                if (!use.safe(iuv, ivj)) continue;

                // Skip SEDs (if asked)
                if ((iuv+ivj) % seds_step != 0) {
                    use.safe(iuv, ivj) = false;
                    continue;
                }

                vec1d tlam = lam.safe(iuv, ivj, _);
                vec1d tsed = sed.safe(iuv, ivj, _);

                tsed = lsun2uJy(zf, df, tlam, tsed);
                tlam *= (1.0 + zf);

                flux.safe(iuv, ivj, 0) = ML_cor*sed2flux(
                    selection_filter.lam, selection_filter.res, tlam, tsed
                );

                for (uint_t l : range(filters)) {
                    flux.safe(iuv, ivj, l+1) = ML_cor*sed2flux(
                        filters[l].lam, filters[l].res, tlam, tsed
                    );
                }
            }

            // Select probability threshold below which SEDs are ignored
            const double pthr = 0.1/count(use);

            // Find minimum mass we need to bother with
            const double mmin = log10(min(flim/flux(_,_,0)));
            const uint_t m0 = lower_bound(m, mmin);

            // Integrate over stellar masses
            for (uint_t im : range(m0, m.size())) {
                const double mm = nm.safe[im];

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
                        on_disk_sed_changed(ised_d);

                        // pdisk(d | z,M*) * N(M*,t,z)
                        const double pdisk = nmz*pblue.safe(iduv, idvj);
                        const double pred_bt1 = nmz*pred.safe(iduv, idvj);

                        // Load stuff
                        const vec1d fdisk = mm*flux.safe(iduv,idvj,_);

                        // Optimization: send B/T=0 now, don't need to iterate over SED bulge
                        if (pdisk >= nmz*pthr) {
                            on_generated(im, it, ised_d, use.safe[0],
                                0, pbt.safe[0]*pdisk, fdisk, fdisk*0.0);
                        }
                        // Optimization: send B/T=1 now, don't need to iterate over SED disk
                        // (trick: pretend we're actually iterating over SED bulge now)
                        if (pred_bt1 >= nmz*pthr) {
                            on_bulge_sed_changed(ised_d);
                            on_generated(im, it, use.safe[0], ised_d,
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
                            on_bulge_sed_changed(ised_b);

                            // p(blue) * pdisk(d | z,M*) * N(M*,t,z)
                            const double pb = pdisk*pblue.safe(ibuv,ibvj);
                            const double pr = pdisk*pred.safe(ibuv,ibvj);

                            // Skip improbable SEDs
                            if (pb < pdisk*pthr && pr < pdisk*pthr) continue;

                            // Load stuff
                            const vec1d fbulge = mm*flux.safe(ibuv,ibvj,_);

                            // Integrate over B/T (skipping B/T=0 and 1)
                            for (uint_t ibt : range(1, bt.size()-1)) {
                                const double btn = bt.safe[ibt];
                                const double bti = 1.0 - btn;

                                // pbulge(b | z,M*,B/T)
                                const double pbulge = (btn >= 0.6 ? pr : 0.5*(pr + pb));

                                // Compute this galaxy's number density (per delta z)
                                const double tprob = pbulge*pbt.safe[ibt];

                                // And forward to derived class.
                                on_generated(im, it, ised_d, ised_b,
                                    ibt, tprob, fdisk*bti, fbulge*btn);
                            }
                        }
                    }
                }
            }
        }
    };
}
