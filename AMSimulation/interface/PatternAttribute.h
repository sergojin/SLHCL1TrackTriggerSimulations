#ifndef AMSimulation_PatternAttribute_h_
#define AMSimulation_PatternAttribute_h_

struct PatternAttribute {
    unsigned long n;  // popularity
    float invPt_mean;
    float invPt_variance;
    float cotTheta_mean;
    float cotTheta_variance;
    float phi_mean;
    float phi_variance;
    float z0_mean;
    float z0_variance;
    //float d0_mean;
    //float d0_variance;

    void reset() {
        n                 = 0;
        invPt_mean        = 0.;
        invPt_variance    = 0.;
        cotTheta_mean     = 0.;
        cotTheta_variance = 0.;
        phi_mean          = 0.;
        phi_variance      = 0.;
        z0_mean           = 0.;
        z0_variance       = 0.;
        //d0_mean           = 0.;
        //d0_variance       = 0.;
    }

    void fill(float invPt, float cotTheta, float phi, float z0) {
        ++n;

        auto fill_statistics = [](float x, unsigned long n, float& mean, float& variance) {
            mean += (x - mean)/n;
            if (n > 1)  variance += (x - mean)*(x - mean)/(n-1) - (variance/n);
        };

        fill_statistics(invPt, n, invPt_mean, invPt_variance);
        fill_statistics(cotTheta, n, cotTheta_mean, cotTheta_variance);
        fill_statistics(phi, n, phi_mean, phi_variance);
        fill_statistics(z0, n, z0_mean, z0_variance);
    }
};

#endif
