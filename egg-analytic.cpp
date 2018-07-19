#include <phypp.hpp>

const std::string egg_share_dir = file::directorize("/home/cschreib/code/egg/share");

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

int phypp_main(int argc, char* argv[]) {
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    uint_t nzbin = zb.size()-1;
    std::string suffix = "rebin2";
    // std::string filters_dir = egg_share_dir+"filter-db/";
    std::string filters_dir = "/home/cschreib/code/egg-analytic/";

    // Define instead a fixed magnitude cut
    double maglim = 24.5;
    std::string selection_band = "euclid-vis";

    // Debug
    constexpr const bool save_psed = false;

    read_args(argc, argv, arg_list(suffix, filters_dir, maglim, selection_band));

    // Internal parameters
    const double flim = e10(0.4*(23.9 - maglim));
    const auto cosmo = get_cosmo("std"); // must be specific cosmology for which EGG is calibrated
    const double area = 1.0*sqr(dpi/180.0)/(4*dpi); // fraction of the sky we simulate

    // Read SED library
    vec2b use;
    vec3f sed, lam;
    vec2f bvj, buv;
    fits::read_table(egg_share_dir+"opt_lib_fast_hd.fits", ftable(use, lam, sed, buv, bvj));

    // Read filter library
    std::string filter_db_file = filters_dir+"db.dat";
    auto filter_db = read_filter_db(filter_db_file);

    filter_t selection_filter;
    if (!get_filter(filter_db, selection_band, selection_filter)) {
        return 1;
    }

    // Mass function
    csmf cs_mf;

    // Base arrays
    const vec1d m  = rgen(4.0,  13.0, 1e2); const double dm  = m[1]  - m[0];
    const vec1d nm = e10(m);
    const vec1d a  = rgen(-5.0, 6.0,  5e1);
    const vec1d b  = rgen(-0.1, 0.2,  5e1);
    const vec1d uv = bin_center(buv);
    const vec1d vj = bin_center(bvj);
    const vec1d bt = rgen(0.0,  1.0,  1e2);

    vec2d map_vj(use.dims), map_uv(use.dims);
    for (uint_t i : range(uv)) {
        map_uv(i,_) = uv[i];
    }
    for (uint_t i : range(vj)) {
        map_vj(_,i) = vj[i];
    }

    // Base distributions
    auto gaussian_integrate = [](double x1, double x2, double mu, double sigma) {
        // Sample (incorect!)
        // return exp(-sqr(x-mu)/(2.0*sqr(sigma)))/(sigma*sqrt(2*dpi));

        // Integrate
        return 0.5*(erf((x2 - mu)/(sigma*sqrt(2.0))) - erf((x1 - mu)/(sigma*sqrt(2.0))));
    };
    auto gaussian = [&](const vec1d& x, double mu, double sigma) {
        vec1d ret(x.dims);
        for (uint_t ix : range(x)) {
            const double x1 = (ix == 0          ? x.safe[ix] : 0.5*(x.safe[ix] + x.safe[ix-1]));
            const double x2 = (ix == x.size()-1 ? x.safe[ix] : 0.5*(x.safe[ix] + x.safe[ix+1]));
            ret.safe[ix] = gaussian_integrate(x1, x2, mu, sigma);
        }

        return ret;
    };
    auto lognormal_integrate = [](double x1, double x2, double mu, double sigma) {
        // Sample (incorect!)
        // return exp(-sqr(log(x/mu))/(2.0*sqr(sigma)))/(x*sigma*sqrt(2*dpi));

        // Integrate
        double up = erf(log(x2/mu)/(sigma*sqrt(2.0)));
        double low = (x1 > 0.0 ? erf(log(x1/mu)/(sqrt(2.0)*sigma)) : -1.0);
        return 0.5*(up - low);
    };
    auto lognormal = [&](const vec1d& x, double mu, double sigma) {
        vec1d ret(x.dims);
        for (uint_t ix : range(x)) {
            const double x1 = (ix == 0          ? x.safe[ix] : 0.5*(x.safe[ix] + x.safe[ix-1]));
            const double x2 = (ix == x.size()-1 ? x.safe[ix] : 0.5*(x.safe[ix] + x.safe[ix+1]));
            ret.safe[ix] = lognormal_integrate(x1, x2, mu, sigma);
        }

        return ret;
    };
    auto upper_lognormal = vectorize_lambda([](double xlow, double mu, double sigma) {
        return 0.5 - 0.5*erf(log(xlow/mu)/(sigma*sqrt(2.0)));
    });
    auto upper_gaussian = vectorize_lambda([](double xlow, double mu, double sigma) {
        return 0.5 - 0.5*erf((xlow - mu)/(sigma*sqrt(2.0)));
    });
    auto lower_gaussian = vectorize_lambda([](double xup, double mu, double sigma) {
        return 0.5 + 0.5*erf((xup - mu)/(sigma*sqrt(2.0)));
    });

    // Average PSF metrics
    vec1d e1(nzbin);
    vec1d e2(nzbin);
    vec1d r2(nzbin);
    vec1d q11(nzbin);
    vec1d q12(nzbin);
    vec1d q22(nzbin);

    // Iterate over redshift bins
    for (uint_t iz : range(nzbin)) {
        std::string zdir = "full_z"+to_string(iz)+"/";
        std::string filename = zdir+"psfs-"+suffix+".txt";
        if (!file::exists(filename)) continue;

        // Read PSF library
        vec1s ssed, zz;
        vec1d me1, me2, mr2, mq11, mq12, mq22;
        ascii::read_table(filename, 0, ssed, zz,
            me1, me2, mr2, mq11, mq12, mq22);

        phypp_check(!zz.empty(), "empty PSF file '", filename, "'");

        // Sort by z then SED
        {
            vec1u ids = sort(zz+ssed);
            ssed = ssed[ids];
            zz  = zz[ids];
            me1 = me1[ids];
            me2 = me2[ids];
            mr2 = mr2[ids];
            mq11 = mq11[ids];
            mq12 = mq12[ids];
            mq22 = mq22[ids];
        }

        // Get SED ID
        vec1u sed_vj(ssed.size()), sed_uv(ssed.size());
        for (uint_t i : range(ssed)) {
            vec1s spl = split(ssed.safe[i], "-");
            from_string(spl.safe[0], sed_uv.safe[i]);
            from_string(spl.safe[1], sed_vj.safe[i]);
        }

        // DEBUG: save p(SED)
        vec2d psed;
        if (save_psed) {
            psed.resize(use.dims);
        }

        const vec1s uz = unique_values_sorted(zz);
        vec1f uzf;
        from_string(replace(uz, "p", "."), uzf);

        const uint_t ntz = uz.size();
        uint_t nsed = 0;

        const vec1d udf = lumdist(uzf, cosmo);

        // Build N(z)
        // double z0 = mean(uzf);
        // vec1d nz = sqr(uzf/z0)*exp(-pow(uzf/z0, 1.5));
        // nz /= integrate(uzf, nz);

        // Integrate over redshift
        auto pg = progress_start(ntz);
        vec1d zq11(ntz), zq12(ntz), zq22(ntz), nz(ntz);
        for (uint_t itz : range(ntz)) {
            const double zf = uzf.safe[itz];
            const double dz = (itz == 0 ? uzf.safe[1] - uzf.safe[0] : uzf.safe[itz] - uzf.safe[itz-1]);
            const double df = udf.safe[itz];

            double tq11 = 0, tq12 = 0, tq22 = 0;
            double tprob = 0;

            vec2d tpsed, tpsed_t, tpsed_m, tpsed_d, tpsed_b, tpsed_bt;
            if (save_psed) {
                tpsed.resize(use.dims);
                tpsed_t.resize(use.dims);
                tpsed_m.resize(use.dims);
                tpsed_d.resize(use.dims);
                tpsed_b.resize(use.dims);
                tpsed_bt.resize(use.dims);
            }

            const auto bz = equal_range(zz, uz[itz]);
            if (nsed == 0) {
                nsed = bz[1]-bz[0];
            } else {
                phypp_check(bz[1]-bz[0] == nsed, "incorrect number of SEDs for z", uz[itz],
                    " (got ", bz[1]-bz[0], ", expected ", nsed, ")");
            }

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

            // Pre-compute fluxes in selection band in library
            vec2d flux(use.dims);
            for (uint_t iuv : range(flux.dims[0]))
            for (uint_t ivj : range(flux.dims[1])) {
                if (!use.safe(iuv, ivj)) continue;

                vec1d tlam = lam.safe(iuv, ivj, _);
                vec1d tsed = sed.safe(iuv, ivj, _);
                tsed = lsun2uJy(zf, df, tlam, tsed)*(selection_filter.rlam/tlam); // photon-weighted flux
                tlam *= (1.0 + zf);

                flux.safe(iuv, ivj) = ML_cor*sed2flux(
                    selection_filter.lam, selection_filter.res, tlam, tsed
                );
            }

            // Find minimum mass we need to bother with
            const double mmin = log10(min(flim/flux));
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

                    double ttq11  = 0;
                    double ttq12  = 0;
                    double ttq22  = 0;
                    double ttprob = 0;

                    // Integrate over disk SED
                    for (uint_t ised_d : range(bz[0], bz[1])) {
                        const uint_t iduv = sed_uv.safe[ised_d];
                        const uint_t idvj = sed_vj.safe[ised_d];

                        // pdisk(d | z,M*)
                        const double pdisk = pblue.safe(iduv, idvj);

                        // Load stuff
                        const double fdisk = mm*flux.safe(iduv,idvj);
                        const double dq11 = mq11.safe[ised_d];
                        const double dq12 = mq12.safe[ised_d];
                        const double dq22 = mq22.safe[ised_d];

                        // Integrate over bulge SED
                        for (uint_t ised_b : range(bz[0], bz[1])) {
                            const uint_t ibuv = sed_uv.safe[ised_b];
                            const uint_t ibvj = sed_vj.safe[ised_b];

                            double pb = pblue.safe(ibuv,ibvj);
                            double pr = pred.safe(ibuv,ibvj);

                            // NB: we factor in the disk p(SED) here to save time
                            pb *= pdisk;
                            pr *= pdisk;

                            // Load stuff
                            const double fbulge = mm*flux.safe(ibuv,ibvj);
                            const double bq11 = mq11.safe[ised_b];
                            const double bq12 = mq12.safe[ised_b];
                            const double bq22 = mq22.safe[ised_b];

                            // Integrate over B/T
                            for (uint_t ibt : range(bt)) {
                                const double btn = bt.safe[ibt];
                                const double bti = 1.0 - btn;

                                // pbulge(b | z,M*,B/T)
                                const double pbulge = (btn >= 0.6 ? pr : 0.5*(pr + pb));

                                // Use full flux for detection limit
                                const double fpred = fbulge*btn + fdisk*bti;

                                if (fpred >= flim) {
                                    // We passed the magnitude cut
                                    const double prob = pbulge*pbt.safe[ibt];

                                    // Compute flux-weighted B/T
                                    const double fbtn = fbulge*btn/fpred;
                                    const double fbti = 1.0 - fbtn;

                                    // Add up PSF
                                    ttq11  += (dq11*fbti + bq11*fbtn)*prob;
                                    ttq12  += (dq12*fbti + bq12*fbtn)*prob;
                                    ttq22  += (dq22*fbti + bq22*fbtn)*prob;
                                    ttprob += prob;

                                    if (save_psed) {
                                        // Save weights
                                        tpsed.safe(iduv,idvj) += fbti*prob;
                                        tpsed.safe(ibuv,ibvj) += fbtn*prob;

                                        tpsed_d.safe(iduv,idvj) += prob;
                                        tpsed_b.safe(ibuv,ibvj) += prob;

                                        // Save type
                                        tpsed_t.safe(iduv,idvj) += it*fbti*prob;
                                        tpsed_t.safe(ibuv,ibvj) += it*fbtn*prob;

                                        // Save mass
                                        tpsed_m.safe(iduv,idvj) += m.safe[im]*fbti*prob;
                                        tpsed_m.safe(ibuv,ibvj) += m.safe[im]*fbtn*prob;

                                        // Save bt
                                        tpsed_bt.safe(iduv,idvj) += fbtn*fbti*prob;
                                        tpsed_bt.safe(ibuv,ibvj) += fbtn*fbtn*prob;
                                    }
                                }
                            }
                        }
                    }

                    // Accumulated probability density for this mass bin
                    // This avoids numerical errors
                    tq11  += ttq11*nmz;
                    tq12  += ttq12*nmz;
                    tq22  += ttq22*nmz;
                    tprob += ttprob*nmz;
                }
            }

            // Average quantities
            tq11 /= tprob;
            tq12 /= tprob;
            tq22 /= tprob;

            // Store
            zq11[itz] = tq11;
            zq12[itz] = tq12;
            zq22[itz] = tq22;
            nz[itz]   = tprob;

            if (save_psed) {
                fits::write(zdir+"psed.fits", tpsed);
                fits::write(zdir+"psed-disk.fits", tpsed_d/total(tpsed_d));
                fits::write(zdir+"psed-bulge.fits", tpsed_b/total(tpsed_b));
                fits::write(zdir+"psed-type.fits", tpsed_t/tpsed);
                fits::write(zdir+"psed-mass.fits", tpsed_m/tpsed);
                fits::write(zdir+"psed-bt.fits", tpsed_bt/tpsed);
            }

            progress(pg);
        }

        // Average over N(z)
        double inz = integrate(uzf, nz);
        q11[iz] = integrate(uzf, zq11*nz)/inz;
        q12[iz] = integrate(uzf, zq12*nz)/inz;
        q22[iz] = integrate(uzf, zq22*nz)/inz;

        // Recompute ellipticity from moments
        e1[iz] = (q11[iz] - q22[iz])/(q11[iz] + q22[iz]);
        e2[iz] = (2*q12[iz])/(q11[iz] + q22[iz]);
        r2[iz] = q11[iz] + q22[iz];

        print("z=[",
            format::fixed(format::precision(zb[iz], 3)), ",",
            format::fixed(format::precision(zb[iz+1], 3)),
            "]: e1=", e1[iz], ", e2=", e2[iz], ", r2=", r2[iz],
            ", Q11=", q11[iz], ", Q12=", q12[iz], ", Q22=", q22[iz],
            " (", ntz, " z steps, ", nsed, " SEDs)");

        fits::write_table(zdir+"psf-mean-"+suffix+".fits",
            "z", uzf, "nz", nz, "q11", zq11, "q12", zq12, "q22", zq22,
            "e1", (zq11-zq22)/(zq11+zq22), "e2", 2*zq12/(zq11+zq22), "r2", zq11+zq22
        );
    }

    fits::write_table("psf-mean-"+suffix+".fits", ftable(zb, e1, e2, r2, q11, q12, q22));

    return 0;
}
