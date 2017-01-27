#!/usr/bin/env python

from rootdrawing import *
from parser import *
from TrackParametersToTT import *

#col  = TColor.GetColor("#1f78b4")  # mu0
#fcol = TColor.GetColor("#a6cee3")  # mu0

#col  = TColor.GetColor("#e31a1c")  # nu140
#fcol = TColor.GetColor("#fb9a99")  # nu140

#col  = TColor.GetColor("#6a3d9a")  # tttt140
#fcol = TColor.GetColor("#cab2d6")  # tttt140

col  = TColor.GetColor("#1A1AE3")  # merged
fcol = TColor.GetColor("#9999FB")  # merged

col1 = TColor.GetColor("#BA0000")
col2 = TColor.GetColor("#CC8000")
col3 = TColor.GetColor("#402033")
fcol1 = TColor.GetColor("#E57D7D")
fcol2 = TColor.GetColor("#FFC766")
fcol3 = TColor.GetColor("#6C5563")


# ______________________________________________________________________________
def drawer_book(options):
    histos = {}

    hname = "nroads_per_event"
    nbins, xmin, xmax = 400, 0., 800.*options.xscale
    histos[hname] = TH1F(hname, "; # roads/tower/BX"                , nbins, xmin, xmax)

    hname = "ncombinations_per_event"
    nbins, xmin, xmax = 400, 0., 1600.*options.xscale
    histos[hname] = TH1F(hname, "; # combinations/tower/BX"         , nbins, xmin, xmax)

    hname = "ncombinations_acb_per_event"
    nbins, xmin, xmax = 400, 0., 1600.*options.xscale
    histos[hname] = TH1F(hname, "; # ACB combinations/tower/BX"     , nbins, xmin, xmax)

    hname = "ntracks_per_event"
    nbins, xmin, xmax = 200, 0., 200.*options.xscale
    histos[hname] = TH1F(hname, "; # tracks/tower/BX"               , nbins, xmin, xmax)

    hname = "ngoods_per_event"
    nbins, xmin, xmax = 200, 0., 200.*options.xscale
    histos[hname] = TH1F(hname, "; # good tracks/tower/BX"          , nbins, xmin, xmax)

    hname = "nduplicates_per_event"
    nbins, xmin, xmax = 200, 0., 200.*options.xscale
    histos[hname] = TH1F(hname, "; # duplicate tracks/tower/BX"     , nbins, xmin, xmax)

    hname = "nfakes_per_event"
    nbins, xmin, xmax = 200, 0., 200.*options.xscale
    histos[hname] = TH1F(hname, "; # fake tracks/tower/BX"          , nbins, xmin, xmax)

    hname = "nparts_per_event"
    nbins, xmin, xmax = 200, 0., 200.*options.xscale
    histos[hname] = TH1F(hname, "; # trkParts/tower/BX"             , nbins, xmin, xmax)

    hname = "nfounds_per_event"
    nbins, xmin, xmax = 200, 0., 200.*options.xscale
    histos[hname] = TH1F(hname, "; # found trkParts/tower/BX"       , nbins, xmin, xmax)

    hname = "foundrate_per_event"
    nbins, xmin, xmax = 120, 0., 1.2
    histos[hname] = TH1F(hname, "; found rate"                      , nbins, xmin, xmax)

    hname = "dupfakerate_per_event"
    nbins, xmin, xmax = 120, 0., 1.2
    histos[hname] = TH1F(hname, "; dup+fake rate"                   , nbins, xmin, xmax)

    hname = "fakerate_per_event"
    nbins, xmin, xmax = 120, 0., 1.2
    histos[hname] = TH1F(hname, "; fake rate"                       , nbins, xmin, xmax)

    for c in ["good", "duplicate", "fake"]:
        hname = "pt_%s" % c
        nbins, xmin, xmax = 80, 0., 80.
        histos[hname] = TH1F(hname, "; track p_{T} [GeV]", nbins, xmin, xmax)

        hname = "eta_%s" % c
        nbins, xmin, xmax = 80, 0., 0.8
        histos[hname] = TH1F(hname, "; track #eta", nbins, xmin, xmax)

        hname = "chi2_%s" % c
        nbins, xmin, xmax = 80, 0., 40
        histos[hname] = TH1F(hname, "; track #chi^{2}/ndof", nbins, xmin, xmax)

    # Change binning
    if options.pu == 0:  # single-track events
        histos["nroads_per_event"           ].SetBins(40, 0., 40.)
        histos["ncombinations_per_event"    ].SetBins(40, 0., 40.)
        histos["ncombinations_acb_per_event"].SetBins(40, 0., 40.)
        histos["ntracks_per_event"          ].SetBins(40, 0., 40.)
        histos["ngoods_per_event"           ].SetBins(40, 0., 40.)
        histos["nduplicates_per_event"      ].SetBins(40, 0., 40.)
        histos["nfakes_per_event"           ].SetBins(40, 0., 40.)
        histos["nparts_per_event"           ].SetBins(40, 0., 40.)
        histos["nfounds_per_event"          ].SetBins(40, 0., 40.)

    # Style
    for hname, h in histos.iteritems():
        h.SetLineWidth(2); h.SetMarkerSize(0)
        h.SetLineColor(col); h.SetFillColor(fcol)
    donotdelete.append(histos)
    return histos

def drawer_project(tree, histos, options):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("trkParts_pt"     , 1)
    tree.SetBranchStatus("trkParts_eta"    , 1)
    tree.SetBranchStatus("trkParts_phi"    , 1)
    #tree.SetBranchStatus("trkParts_vx"     , 1)
    #tree.SetBranchStatus("trkParts_vy"     , 1)
    tree.SetBranchStatus("trkParts_vz"     , 1)
    tree.SetBranchStatus("trkParts_charge" , 1)
    tree.SetBranchStatus("trkParts_primary", 1)
    tree.SetBranchStatus("trkParts_intime" , 1)
    tree.SetBranchStatus("trkParts_signal" , 1)
    tree.SetBranchStatus("trkParts_pdgId"  , 1)
    tree.SetBranchStatus("AMTTRoads_patternRef"  , 1)
    tree.SetBranchStatus("AMTTRoads_stubRefs"    , 1)
    #tree.SetBranchStatus("AMTTTracks_invPt"      , 1)
    #tree.SetBranchStatus("AMTTTracks_phi0"       , 1)
    #tree.SetBranchStatus("AMTTTracks_cottheta"   , 1)
    #tree.SetBranchStatus("AMTTTracks_z0"         , 1)
    tree.SetBranchStatus("AMTTTracks_pt"         , 1)
    tree.SetBranchStatus("AMTTTracks_eta"        , 1)
    tree.SetBranchStatus("AMTTTracks_chi2"       , 1)
    tree.SetBranchStatus("AMTTTracks_ndof"       , 1)
    tree.SetBranchStatus("AMTTTracks_synMatchCat", 1)
    tree.SetBranchStatus("AMTTTracks_synTpId"    , 1)
    tree.SetBranchStatus("AMTTTracks_patternRef" , 1)
    tree.SetBranchStatus("AMTTTracks_roadRef"    , 1)
    #tree.SetBranchStatus("AMTTTracks_stubRefs"   , 1)

    # Loop over events
    for ievt, evt in enumerate(tree):
        if (ievt == options.nentries):  break

        if (ievt % 100 == 0):  print "Processing event: %i" % ievt

        # Loop over tracking particles
        nparts_all = evt.trkParts_primary.size()
        trkparts = {}
        if options.no_tp:  nparts_all = 0

        for ipart in xrange(nparts_all):
            if options.pu == 0:  # single-track events
                if ipart > 0:
                    break

            charge  = evt.trkParts_charge [ipart]
            primary = evt.trkParts_primary[ipart]
            intime  = evt.trkParts_intime [ipart]
            signal  = evt.trkParts_signal [ipart]
            pt      = evt.trkParts_pt     [ipart]

            if not (charge!=0 and primary and intime and pt >= 1.):
                continue

            if options.signal and not signal:
                continue

            eta     = evt.trkParts_eta    [ipart]
            phi     = evt.trkParts_phi    [ipart]
            #vx      = evt.trkParts_vx     [ipart]
            #vy      = evt.trkParts_vy     [ipart]
            vz      = evt.trkParts_vz     [ipart]
            pdgId   = evt.trkParts_pdgId  [ipart]

            aux_TT = TrackParametersToTT(phi, float(charge)/pt , eta, vz, apply_pt_cut=False)
            if aux_TT != options.tower:
                continue

            trkparts[ipart] = 0
            if options.verbose:  print ievt, "part ", ipart, pt, eta, phi, vz

        # Loop over reconstructed roads
        nroads_all = evt.AMTTRoads_patternRef.size()
        nroads, ncombs, ncombs_acb = 0, 0, 0

        for iroad in xrange(nroads_all):
            patternRef = evt.AMTTRoads_patternRef[iroad]
            if not (patternRef < options.maxPatterns):
                continue

            roadRef = iroad
            if not (roadRef < options.maxRoads):
                continue

            nroads += 1

            stubRefs   = evt.AMTTRoads_stubRefs  [iroad]
            assert(stubRefs.size() == 6)

            ncombs_per_road = 1
            for v in stubRefs:
                n = v.size()
                if n != 0:
                    ncombs_per_road *= n
            ncombs += ncombs_per_road

            ncombs_acb_per_road = ncombs_per_road
            is_6oo6 = (sum([v.size() > 0 for v in stubRefs]) == 6)
            if is_6oo6:
                ncombs_acb_per_road += (stubRefs[0].size() * stubRefs[1].size() * stubRefs[2].size() * stubRefs[3].size() * stubRefs[4].size())
                ncombs_acb_per_road += (stubRefs[0].size() * stubRefs[1].size() * stubRefs[2].size() * stubRefs[3].size() * stubRefs[5].size())
                ncombs_acb_per_road += (stubRefs[0].size() * stubRefs[1].size() * stubRefs[2].size() * stubRefs[4].size() * stubRefs[5].size())
                ncombs_acb_per_road += (stubRefs[0].size() * stubRefs[1].size() * stubRefs[3].size() * stubRefs[4].size() * stubRefs[5].size())
                ncombs_acb_per_road += (stubRefs[0].size() * stubRefs[2].size() * stubRefs[3].size() * stubRefs[4].size() * stubRefs[5].size())
                ncombs_acb_per_road += (stubRefs[1].size() * stubRefs[2].size() * stubRefs[3].size() * stubRefs[4].size() * stubRefs[5].size())
            ncombs_acb += ncombs_acb_per_road

        # Loop over reconstructed tracks
        ntracks_all = evt.AMTTTracks_patternRef.size()
        ntracks, ngoods, nduplicates, nfakes = 0, 0, 0, 0

        for itrack in xrange(ntracks_all):
            patternRef  = evt.AMTTTracks_patternRef [itrack]
            if not (patternRef < options.maxPatterns):
                continue

            roadRef     = evt.AMTTTracks_roadRef    [itrack]
            if not (roadRef < options.maxRoads):
                continue

            synMatchCat = evt.AMTTTracks_synMatchCat[itrack]
            synTpId     = evt.AMTTTracks_synTpId    [itrack]
            track_pt    = evt.AMTTTracks_pt         [itrack]
            track_eta   = evt.AMTTTracks_eta        [itrack]
            track_chi2  = evt.AMTTTracks_chi2       [itrack]
            track_ndof  = evt.AMTTTracks_ndof       [itrack]

            ntracks += 1

            if synMatchCat == -2:
                histos["pt_fake"  ].Fill(track_pt)
                histos["eta_fake" ].Fill(track_eta)
                histos["chi2_fake"].Fill(track_chi2/track_ndof)
                nfakes += 1

            elif synMatchCat == -1:
                histos["pt_duplicate" ].Fill(track_pt)
                histos["eta_duplicate"].Fill(track_eta)
                histos["chi2_duplicate"].Fill(track_chi2/track_ndof)
                nduplicates += 1

            elif synMatchCat == 1:
                histos["pt_good"  ].Fill(track_pt)
                histos["eta_good" ].Fill(track_eta)
                histos["chi2_good"].Fill(track_chi2/track_ndof)
                ngoods += 1
                if synTpId in trkparts:
                    trkparts[synTpId] = 1

            else:
                raise Exception("Unexpected synMatchCat: %i" % synMatchCat)

        nparts = len(trkparts)
        nfounds = sum([v for k, v in trkparts.iteritems()])

        if options.verbose:  print ievt, nroads, ncombs, ntracks, ngoods, nduplicates, nfakes, nparts, nfounds

        assert(ntracks == ngoods + nduplicates + nfakes)
        histos["nroads_per_event"           ].Fill(nroads)
        histos["ncombinations_per_event"    ].Fill(ncombs)
        histos["ncombinations_acb_per_event"].Fill(ncombs_acb)
        histos["ntracks_per_event"          ].Fill(ntracks)
        histos["ngoods_per_event"           ].Fill(ngoods)
        histos["nduplicates_per_event"      ].Fill(nduplicates)
        histos["nfakes_per_event"           ].Fill(nfakes)

        assert(nfounds <= nparts)
        histos["nparts_per_event"           ].Fill(nparts)
        histos["nfounds_per_event"          ].Fill(nfounds)

        if nparts:
            histos["foundrate_per_event"        ].Fill(float(nfounds)/nparts)
        if ntracks:
            histos["dupfakerate_per_event"      ].Fill(float(nduplicates+nfakes)/(nduplicates+nfakes+ngoods))
            histos["fakerate_per_event"         ].Fill(float(nfakes)/(nfakes+ngoods))

    tree.SetBranchStatus("*", 1)
    return

def drawer_draw(histos, options):
    options.logy = True

    def displayQuantiles(h, in_quantiles=[0.95,0.99], scalebox=(1.,1.)):
        # Display one-sided confidence intervals, a.k.a quantiles
        n = len(in_quantiles)
        in_quantiles = array('d', in_quantiles)
        quantiles = array('d', [0. for i in xrange(n)])
        h.GetQuantiles(n, quantiles, in_quantiles)

        gPad.Modified(); gPad.Update()
        ps = h.FindObject("stats")
        ps.SetName("mystats")

        newX1NDC = ps.GetX2NDC() - (ps.GetX2NDC() - ps.GetX1NDC()) * scalebox[0]
        newY1NDC = ps.GetY2NDC() - ((ps.GetY2NDC() - ps.GetY1NDC()) / 5 * (5 + n)) * scalebox[1]
        ps.SetX1NDC(newX1NDC)
        ps.SetY1NDC(newY1NDC)

        for iq, q in enumerate(in_quantiles):
            ps.AddText("%i%% CI = %6.4g" % (int(q*100), quantiles[iq]))
        h.stats = [h.GetMean()] + quantiles.tolist()

        h.SetStats(0)
        #gPad.Modified(); gPad.Update()
        ps.Draw()

    for hname, h in histos.iteritems():
        if h.ClassName() == "TH1F":
            if options.logy:
                h.SetMaximum(h.GetMaximum() * 14); h.SetMinimum(0.5)
            else:
                h.SetMaximum(h.GetMaximum() * 1.4); h.SetMinimum(0.)

            h.Draw("hist")
            gPad.SetLogy(options.logy)
            if hname.endswith("_per_event"):
                displayQuantiles(h)

        CMS_label()
        save(options.outdir, "%s_%s" % (hname, options.ss), dot_root=True)
    return

def drawer_draw2(histos, options):
    options.logy = False

    # Specialized
    for v in ["pt", "eta", "chi2"]:
        hname1 = "%s_good" % v
        hname2 = "%s_duplicate" % v
        hname3 = "%s_fake" % v

        h1 = histos[hname1]
        h2 = histos[hname2]
        h3 = histos[hname3]

        h1.SetLineColor(col1); h1.SetMarkerColor(col1); h1.SetFillColor(fcol1)
        h2.SetLineColor(col2); h2.SetMarkerColor(col2); h2.SetFillColor(fcol2)
        h3.SetLineColor(col3); h3.SetMarkerColor(col3); h3.SetFillColor(fcol3)

        # Stack
        hstack1 = h1.Clone(hname1 + "_stack")
        hstack2 = h2.Clone(hname2 + "_stack")
        hstack3 = h3.Clone(hname3 + "_stack")

        hstack2.Add(hstack1)
        hstack3.Add(hstack2)

        hstack3.SetMaximum(hstack3.GetMaximum()*1.4); hstack3.SetMinimum(0)
        hstack3.SetStats(0); hstack3.Draw("hist")
        hstack2.Draw("hist same")
        hstack1.Draw("hist same")
        gPad.SetLogy(options.logy)

        moveLegend(0.64,0.82,0.94,0.94); tlegend.Clear()
        tlegend.AddEntry(h1, "good", "f")
        tlegend.AddEntry(h2, "duplicate", "f")
        tlegend.AddEntry(h3, "fake", "f")
        tlegend.Draw()

        CMS_label()
        save(options.outdir, "%s_stack_%s" % (v, options.ss), dot_root=True)

        # Ratio
        hratio1 = hstack1.Clone(hname1 + "_ratio")
        hratio2 = hstack2.Clone(hname2 + "_ratio")
        hratio3 = hstack3.Clone(hname3 + "_ratio")

        hratio1.Divide(hstack3)
        hratio2.Divide(hstack3)
        hratio3.Divide(hstack3)

        hratio3.SetMaximum(1.2); hratio3.SetMinimum(0)
        hratio3.SetStats(0); hratio3.Draw("hist")
        hratio2.SetMaximum(1.2); hratio2.SetMinimum(0)
        hratio2.Draw("hist same")
        hratio1.SetMaximum(1.2); hratio1.SetMinimum(0)
        hratio1.Draw("hist same")
        gPad.SetLogy(options.logy)

        moveLegend(0.64,0.82,0.94,0.94); tlegend.Clear()
        tlegend.AddEntry(h1, "good", "f")
        tlegend.AddEntry(h2, "duplicate", "f")
        tlegend.AddEntry(h3, "fake", "f")
        tlegend.Draw()

        CMS_label()
        save(options.outdir, "%s_ratio_%s" % (v, options.ss), dot_root=True)

        # Norm
        hnorm1 = h1.Clone(hname1 + "_norm")
        hnorm2 = h2.Clone(hname2 + "_norm")
        hnorm3 = h3.Clone(hname3 + "_norm")

        hnorm1.Scale(1.0/hnorm1.Integral())
        hnorm2.Scale(1.0/hnorm2.Integral())
        hnorm3.Scale(1.0/hnorm3.Integral())

        hnorm1.SetFillStyle(0)
        hnorm2.SetFillStyle(0)
        hnorm3.SetFillStyle(0)

        ymax = max(hnorm1.GetMaximum(), hnorm2.GetMaximum(), hnorm3.GetMaximum())
        hnorm3.SetMaximum(ymax*1.4); hnorm3.SetMinimum(0)
        hnorm3.SetStats(0); hnorm3.Draw("hist")
        hnorm2.Draw("hist same")
        hnorm1.Draw("hist same")
        gPad.SetLogy(options.logy)

        moveLegend(0.64,0.82,0.94,0.94); tlegend.Clear()
        tlegend.AddEntry(h1, "good", "l")
        tlegend.AddEntry(h2, "duplicate", "l")
        tlegend.AddEntry(h3, "fake", "l")
        tlegend.Draw()

        CMS_label()
        save(options.outdir, "%s_norm_%s" % (v, options.ss), dot_root=True)

        # Purity
        hpurity1 = h1.Clone(hname1 + "_purity")
        hpurity2 = h2.Clone(hname2 + "_purity")
        hpurity3 = h3.Clone(hname3 + "_purity")

        hpurity1 = hpurity1     # num = good
        hpurity3.Add(hpurity1)  # denom = good + fake

        hpurity1.Divide(hpurity3)
        hpurity1.SetMaximum(1.2); hpurity1.SetMinimum(0)
        hpurity1.SetStats(0); hpurity1.Draw("p")

        CMS_label()
        save(options.outdir, "%s_purity_%s" % (v, options.ss), dot_root=True)

        donotdelete.append([hstack1, hstack2, hstack3])
        donotdelete.append([hratio1, hratio2, hratio3])
        donotdelete.append([hnorm1, hnorm2, hnorm3])
        donotdelete.append([hpurity1, hpurity2, hpurity3])
    return

def drawer_sitrep(histos, options):
    print "--- SITREP ---------------------------------------------------------"
    print "--- Using tt{0}, pu{1}, ss={2}, maxPatterns={3}".format(options.tower, options.pu, options.ss, options.maxPatterns)
    print "--- Variable, mean, 95%% CI, 99%% CI:"
    h = histos["nroads_per_event"]
    print "nroads per event\n{0:6.4g}\n{1:6.4g}\n{2:6.4g}\n".format(*h.stats)
    h = histos["ncombinations_per_event"]
    print "ncombs per event\n{0:6.4g}\n{1:6.4g}\n{2:6.4g}\n".format(*h.stats)
    h = histos["ntracks_per_event"]
    print "ntrks  per event\n{0:6.4g}\n{1:6.4g}\n{2:6.4g}\n".format(*h.stats)
    #h = histos["ngoods_per_event"]
    #print "ngoods per event\n{0:6.4g}\n{1:6.4g}\n{2:6.4g}".format(*h.stats)
    #h = histos["nduplicates_per_event"]
    #print "ndupls per event\n{0:6.4g}\n{1:6.4g}\n{2:6.4g}".format(*h.stats)
    #h = histos["nfakes_per_event"]
    #print "nfakes per event\n{0:6.4g}\n{1:6.4g}\n{2:6.4g}".format(*h.stats)

    print "--- Efficiency, dup+fake rate, fake rate:"
    foundrate   = histos["foundrate_per_event"].GetMean()
    dupfakerate = histos["dupfakerate_per_event"].GetMean()
    fakerate    = histos["fakerate_per_event"].GetMean()
    print "eff vs mis-id\n{0:6.4g}\n{1:6.4g}\n{2:6.4g}\n".format(foundrate, dupfakerate, fakerate)


# ______________________________________________________________________________
# Main function
def main(options):

    # Init
    drawerInit = DrawerInit()
    tchain = TChain("ntupler/tree", "")
    tchain.AddFileInfoList(options.tfilecoll.GetList())

    # Process
    histos = drawer_book(options)
    drawer_project(tchain, histos, options)
    drawer_draw(histos, options)
    drawer_draw2(histos, options)
    drawer_sitrep(histos, options)


# ______________________________________________________________________________
if __name__ == '__main__':

    # Setup argument parser
    parser = argparse.ArgumentParser()

    # Add default arguments
    add_drawer_arguments(parser)

    # Add more arguments
    parser.add_argument("ss", help="short name of superstrip definition (e.g. ss256)")
    parser.add_argument("--maxPatterns", type=int, default=999999999, help="number of patterns to reach the desired coverage")
    parser.add_argument("--maxRoads", type=int, default=999999999, help="number of roads to not truncate")
    parser.add_argument("--minPt", type=float, default=3, help="min pT for gen particle (default: %(default)s)")
    parser.add_argument("--no-tp", action="store_true", help="ignore tracking particles (default: %(default)s)")
    parser.add_argument("--xscale", type=float, default=1, help="scale factor for the x-axis range (default: %(default)s)")
    parser.set_defaults(pu=140)

    # Parse default arguments
    options = parser.parse_args()
    parse_drawer_options(options)
    options.ptmin = options.minPt

    # Call the main function
    main(options)
