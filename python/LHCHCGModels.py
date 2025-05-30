import ROOT
from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder

## Naming conventions
CMS_to_LHCHCG_Dec = {
    "hww": "WW",
    "hzz": "ZZ",
    "hgg": "gamgam",
    "hbb": "bb",
    "hcc": "cc",
    "htt": "tautau",
    "hmm": "mumu",
    "hzg": "Zgam",
    "hgluglu": "gluglu",
    "hinv": "inv",
}
CMS_to_LHCHCG_DecSimple = {
    "hww": "WW",
    "hzz": "ZZ",
    "hgg": "gamgam",
    "hbb": "bb",
    "hcc": "bb",
    "htt": "tautau",
    "hmm": "mumu",
    "hzg": "gamgam",
    "hgluglu": "bb",
    "hinv": "inv",
}
CMS_to_LHCHCG_Prod = {
    "ggH": "ggF",
    "qqH": "VBF",
    "WH": "WH",
    "WPlusH": "WH",
    "WMinusH": "WH",
    "ZH": "qqZH",
    "ggZH": "ggZH",
    "ttH": "ttH",
    "tHq": "tHjb",
    "tHW": "WtH",
    "bbH": "bbH",
}


class LHCHCGBaseModel(SMLikeHiggsModel):
    def __init__(self):
        SMLikeHiggsModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False
        self.add_bbH = []
        self.bbH_pdf = "pdf_Higgs_gg"
        self.promote_hmm = True
        self.promote_hccgluglu = False
        self.promote_hzg = False

    def preProcessNuisances(self, nuisances):
        if self.add_bbH and not any(row for row in nuisances if row[0] == "QCDscale_bbH"):
            nuisances.append(("QCDscale_bbH", False, "param", ["0", "1"], []))
        if self.add_bbH and not any(row for row in nuisances if row[0] == self.bbH_pdf):
            nuisances.append((bbH_pdf, False, "param", ["0", "1"], []))

    def setPhysicsOptionsBase(self, physOptions):
        for po in physOptions:
            if po.startswith("dohmm="):
                self.promote_hmm = po.replace("dohmm=", "") in [
                    "yes",
                    "1",
                    "Yes",
                    "True",
                    "true",
                ]
                if not self.promote_hmm:
                    print("WARNING: hmm signal will be treated as htt")
                    CMS_to_LHCHCG_DecSimple["hmm"] = "tautau"
            if po.startswith("dohzg="):
                self.promote_hzg = po.replace("dohzg=", "") in [
                    "yes",
                    "1",
                    "Yes",
                    "True",
                    "true",
                ]
                if self.promote_hzg:
                    print("WARNING: hzg signal will be treated standalone")
                    CMS_to_LHCHCG_DecSimple["hzg"] = "Zgam"
            if po.startswith("bbh="):
                self.add_bbH = [d.strip() for d in po.replace("bbh=", "").split(",")]
            if po.startswith("dohccgluglu="):
                self.promote_hccgluglu = po.replace("dohccgluglu=", "") in [
                    "yes",
                    "1",
                    "Yes",
                    "True",
                    "true",
                ]
                if self.promote_hccgluglu:
                    print("Treating hcc and hgg as independent processes")
                    CMS_to_LHCHCG_DecSimple["hcc"] = "cc"
                    CMS_to_LHCHCG_DecSimple["hgluglu"] = "gluglu"
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=", "").split(",")
                print("The Higgs mass range:", self.mHRange)
                if len(self.mHRange) != 2:
                    raise RuntimeError("Higgs mass range definition requires two extrema")
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError("Extrema for Higgs mass range defined with inverterd order. Second must be larger than the first")
        print("Will add bbH to signals in the following Higgs boson decay modes: %s" % (", ".join(self.add_bbH)))

    def dobbH(self):
        self.modelBuilder.doVar("QCDscale_bbH[-7,7]")
        self.modelBuilder.doVar(self.bbH_pdf + "[-7,7]")
        scaler7 = ROOT.ProcessNormalization("CMS_bbH_scaler_7TeV", "", 1.0)
        scaler8 = ROOT.ProcessNormalization("CMS_bbH_scaler_8TeV", "", 1.0)
        scaler13 = ROOT.ProcessNormalization("CMS_bbH_scaler_13TeV", "", 1.0)
        self.modelBuilder.out.safe_import(scaler7)
        self.modelBuilder.out.safe_import(scaler8)
        self.modelBuilder.out.safe_import(scaler13)
        ## YR3
        # self.modelBuilder.out.function("CMS_bbH_scaler_7TeV").addAsymmLogNormal(1.0/1.145, 1.106, self.modelBuilder.out.var("QCDscale_bbH"))
        # self.modelBuilder.out.function("CMS_bbH_scaler_8TeV").addAsymmLogNormal(1.0/1.148, 1.103, self.modelBuilder.out.var("QCDscale_bbH"))
        # self.modelBuilder.out.function("CMS_bbH_scaler_13TeV").addAsymmLogNormal(1.0/1.148, 1.103, self.modelBuilder.out.var("QCDscale_bbH"))
        ## YR4: QCDscale+pdf+alphas
        self.modelBuilder.out.function("CMS_bbH_scaler_7TeV").addAsymmLogNormal(1.0 / 1.225, 1.207, self.modelBuilder.out.var("QCDscale_bbH"))
        self.modelBuilder.out.function("CMS_bbH_scaler_8TeV").addAsymmLogNormal(1.0 / 1.224, 1.206, self.modelBuilder.out.var("QCDscale_bbH"))
        self.modelBuilder.out.function("CMS_bbH_scaler_13TeV").addAsymmLogNormal(1.0 / 1.239, 1.205, self.modelBuilder.out.var("QCDscale_bbH"))
        ## pdf split is reported in YR3. pdf+alphas. I didn't subtract them from above.
        self.modelBuilder.out.function("CMS_bbH_scaler_7TeV").addLogNormal(1.061, self.modelBuilder.out.var(self.bbH_pdf))
        self.modelBuilder.out.function("CMS_bbH_scaler_8TeV").addLogNormal(1.062, self.modelBuilder.out.var(self.bbH_pdf))
        self.modelBuilder.out.function("CMS_bbH_scaler_13TeV").addLogNormal(1.061, self.modelBuilder.out.var(self.bbH_pdf))

    def doMH(self):
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]), float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar(f"MH[{self.mHRange[0]},{self.mHRange[1]}]")
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)


class SignalStrengths(LHCHCGBaseModel):
    "Allow different signal strength fits"

    def __init__(self):
        LHCHCGBaseModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it
        self.POIs = "mu"

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("poi="):
                self.POIs = po.replace("poi=", "")

    def doVar(self, x, constant=True):
        self.modelBuilder.doVar(x)
        vname = re.sub(r"\[.*", "", x)
        self.modelBuilder.out.var(vname).setConstant(constant)
        print(f"SignalStrengths:: declaring {vname} as {x}, and set to constant")

    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.doMH()
        self.doVar("mu[1,0,5]")
        for X in CMS_to_LHCHCG_Dec.values():
            self.doVar("mu_BR_%s[1,0,5]" % X)
            self.doVar("mu_BR7_%s[1,0,5]" % X)
            self.doVar("mu_BR8_%s[1,0,5]" % X)
            self.doVar("mu_BR13_%s[1,0,5]" % X)
        for X in list(CMS_to_LHCHCG_Prod.values()) + [
            "ZH",
            "tH",
            "ggFbbH",
            "ttHtH",
            "VH",
            "VHttHtH",
        ]:
            self.doVar("mu_XS_%s[1,0,5]" % X)
            self.doVar("mu_XS7_%s[1,0,5]" % X)
            self.doVar("mu_XS8_%s[1,0,5]" % X)
            self.doVar("mu_XS13_%s[1,0,5]" % X)
            for Y in CMS_to_LHCHCG_Dec.values():
                self.doVar(f"mu_XS_{X}_BR_{Y}[1,0,5]")
        for X in CMS_to_LHCHCG_DecSimple.values():
            self.doVar("mu_V_%s[1,0,5]" % X)
            self.doVar("mu_F_%s[1,0,5]" % X)
        print("Default parameters of interest: ", self.POIs)
        self.modelBuilder.doSet("POI", self.POIs)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        for P in ALL_HIGGS_PROD:
            if P == "VH":
                continue  # skip aggregates
            for D in SM_HIGG_DECAYS:
                for E in 7, 8, 13:
                    DS = CMS_to_LHCHCG_DecSimple[D]
                    terms = [
                        "mu",
                        "mu_BR_" + CMS_to_LHCHCG_DecSimple[D],
                        "mu_BR" + str(E) + "_" + CMS_to_LHCHCG_DecSimple[D],
                    ]
                    # Hack for ggH
                    if D in self.add_bbH and P == "ggH":
                        b2g = "CMS_R_bbH_ggH_%s_%dTeV[%g]" % (D, E, 0.01)
                        ggs = ",".join(["mu_XS_ggF", "mu_XS%d_ggF" % E])
                        bbs = ",".join(["mu_XS_bbH", "mu_XS%d_bbH" % E, "CMS_bbH_scaler_%dTeV" % E])
                        ## FIXME should include the here also logNormal for QCDscale_bbH
                        self.modelBuilder.factory_('expr::ggH_bbH_sum_%s_%dTeV("@1*@2+@0*@3*@4*@5",%s,%s,%s)' % (D, E, b2g, ggs, bbs))
                        terms += [
                            "ggH_bbH_sum_%s_%dTeV" % (D, E),
                            "mu_XS_ggFbbH",
                            "mu_XS%d_ggFbbH" % E,
                        ]
                        terms += ["mu_XS_ggFbbH_BR_%s" % DS]
                    else:
                        if P in ["ggH", "bbH"]:
                            terms += ["mu_XS_ggFbbH", "mu_XS%d_ggFbbH" % E]
                            terms += ["mu_XS_ggFbbH_BR_%s" % DS]
                        terms += [
                            "mu_XS_" + CMS_to_LHCHCG_Prod[P],
                            "mu_XS%d_%s" % (E, CMS_to_LHCHCG_Prod[P]),
                        ]
                        terms += ["mu_XS_" + CMS_to_LHCHCG_Prod[P] + "_BR_%s" % DS]

                    # Summary modes
                    if P in ["tHW", "tHq"]:
                        terms += ["mu_XS_tH", "mu_XS%d_tH" % E]
                        terms += ["mu_XS_tH_BR_%s" % DS]
                    if P in ["tHW", "tHq", "ttH"]:
                        terms += ["mu_XS_ttHtH", "mu_XS%d_ttHtH" % E]
                        terms += ["mu_XS_ttHtH_BR_%s" % DS]
                    if P in ["ggZH", "ZH"]:
                        terms += ["mu_XS_ZH", "mu_XS%d_ZH" % E]
                        terms += ["mu_XS_ZH_BR_%s" % DS]
                    if P in ["WH", "ZH", "ggZH"]:
                        terms += ["mu_XS_VH", "mu_XS%d_VH" % E]
                        terms += ["mu_XS_VH_BR_%s" % DS]
                    if P in ["WH", "ZH", "ggZH", "tHW", "tHq", "ttH"]:
                        terms += ["mu_XS_VHttHtH", "mu_XS%d_VHttHtH" % E]
                        terms += ["mu_XS_VHttHtH_BR_%s" % DS]
                    # for 2D scans
                    if P in ["ggH", "ttH", "bbH", "tHq", "tHW"]:
                        terms += ["mu_F_" + CMS_to_LHCHCG_DecSimple[D]]
                    else:
                        terms += ["mu_V_" + CMS_to_LHCHCG_DecSimple[D]]
                    self.modelBuilder.factory_("prod::scaling_%s_%s_%dTeV(%s)" % (P, D, E, ",".join(terms)))
                    self.modelBuilder.out.function("scaling_%s_%s_%dTeV" % (P, D, E)).Print("")

    def getHiggsSignalYieldScale(self, production, decay, energy):
        return f"scaling_{production}_{decay}_{energy}"


class SignalStrengthRatios(LHCHCGBaseModel):
    "Allow for fits of ratios of signal strengths"

    def __init__(self):
        LHCHCGBaseModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("poi="):
                self.POIs = po.replace("poi=", "")

    def doVar(self, x, constant=True):
        self.modelBuilder.doVar(x)
        vname = re.sub(r"\[.*", "", x)
        self.modelBuilder.out.var(vname).setConstant(constant)
        print(f"SignalStrengthRatios:: declaring {vname} as {x}, and set to constant")

    def doParametersOfInterest(self):
        """Create POI out of signal strength ratios and MH"""
        self.doMH()
        self.doVar("mu_V_r_F[1,0,5]")
        for X in CMS_to_LHCHCG_DecSimple.values():
            self.doVar("mu_F_%s[1,0,5]" % X)
            self.doVar("mu_V_r_F_%s[1,0,5]" % X)
        self.POIs = "mu_V_r_F"
        print("Default parameters of interest: ", self.POIs)
        self.modelBuilder.doSet("POI", self.POIs)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        for P in ALL_HIGGS_PROD:
            if P == "VH":
                continue  # skip aggregates
            for D in SM_HIGG_DECAYS:
                for E in 7, 8, 13:
                    terms = ["mu_F_" + CMS_to_LHCHCG_DecSimple[D]]
                    if P in ["qqH", "VBF", "VH", "WH", "ZH", "qqZH", "ggZH"]:
                        terms += ["mu_V_r_F", "mu_V_r_F_" + CMS_to_LHCHCG_DecSimple[D]]

                    # Hack for ggH
                    if D in self.add_bbH and P == "ggH":
                        b2g = "CMS_R_bbH_ggH_%s_%dTeV[%g]" % (D, E, 0.01)
                        bbs = "CMS_bbH_scaler_%dTeV" % E
                        self.modelBuilder.factory_('expr::ggH_bbH_sum_%s_%dTeV("1+@0*@1",%s,%s)' % (D, E, b2g, bbs))
                        terms += ["ggH_bbH_sum_%s_%dTeV" % (D, E)]

                    self.modelBuilder.factory_("prod::scaling_%s_%s_%dTeV(%s)" % (P, D, E, ",".join(terms)))
                    self.modelBuilder.out.function("scaling_%s_%s_%dTeV" % (P, D, E)).Print("")

    def getHiggsSignalYieldScale(self, production, decay, energy):
        return f"scaling_{production}_{decay}_{energy}"


class XSBRratios(LHCHCGBaseModel):
    "Model with the ratios of cross sections and branching ratios"

    def __init__(self, denominator="WW"):
        LHCHCGBaseModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it
        self.denominator = denominator

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("poi="):
                self.POIs = po.replace("poi=", "")

    def doVar(self, x, constant=True):
        self.modelBuilder.doVar(x)
        vname = re.sub(r"\[.*", "", x)
        self.modelBuilder.out.var(vname).setConstant(constant)
        print(f"XSBRratios:: declaring {vname} as {x}, and set to constant")

    def doParametersOfInterest(self):
        """Create POI out of signal strength ratios and MH"""
        self.doMH()
        self.doVar("mu_XS7_r_XS8_ggF[1,0,5]")
        self.doVar("mu_XS7_r_XS8_VBF[1,0,5]")
        self.doVar("mu_XS7_r_XS8_WH[1,0,5]")
        self.doVar("mu_XS7_r_XS8_ZH[1,0,5]")
        self.doVar("mu_XS7_r_XS8_ttH[1,0,5]")
        self.modelBuilder.doVar("mu_XS_ggF_x_BR_%s[1,0,5]" % self.denominator)
        self.modelBuilder.doVar("mu_XS_VBF_r_XS_ggF[1,0,5]")
        self.modelBuilder.doVar("mu_XS_WH_r_XS_ggF[1,0,5]")
        self.modelBuilder.doVar("mu_XS_ZH_r_XS_ggF[1,0,5]")
        self.modelBuilder.doVar("mu_XS_ttH_r_XS_ggF[1,0,5]")
        self.POIs = "mu_XS_ggF_x_BR_%s,mu_XS_VBF_r_XS_ggF,mu_XS_ttH_r_XS_ggF,mu_XS_WH_r_XS_ggF,mu_XS_ZH_r_XS_ggF" % self.denominator
        for X in ["ZZ", "tautau", "bb", "gamgam", "WW", "mumu"]:
            if X == self.denominator:
                continue
            self.modelBuilder.doVar(f"mu_BR_{X}_r_BR_{self.denominator}[1,0,5]")
            self.POIs += f",mu_BR_{X}_r_BR_{self.denominator}"
        print("Default parameters of interest: ", self.POIs)
        self.modelBuilder.doSet("POI", self.POIs)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        for P in ALL_HIGGS_PROD:
            if P == "VH":
                continue  # skip aggregates
            for D in SM_HIGG_DECAYS:
                for E in 7, 8, 13:
                    terms = ["mu_XS_ggF_x_BR_%s" % self.denominator]
                    if CMS_to_LHCHCG_DecSimple[D] != self.denominator:
                        terms += [f"mu_BR_{CMS_to_LHCHCG_DecSimple[D]}_r_BR_{self.denominator}"]
                    if P in ["ttH", "tHq", "tHW"]:
                        terms += ["mu_XS_ttH_r_XS_ggF"]
                    if P in ["ZH", "ggZH"]:
                        terms += ["mu_XS_ZH_r_XS_ggF"]
                    if P in ["WH"]:
                        terms += ["mu_XS_WH_r_XS_ggF"]
                    if P in ["qqH"]:
                        terms += ["mu_XS_VBF_r_XS_ggF"]
                    if E == 7:
                        if P in ["ggH", "bbH"]:
                            terms += ["mu_XS7_r_XS8_ggF"]
                        if P in ["ttH", "tHq", "tHW"]:
                            terms += ["mu_XS7_r_XS8_ttH"]
                        if P in ["ZH", "ggZH"]:
                            terms += ["mu_XS7_r_XS8_ZH"]
                        if P in ["WH"]:
                            terms += ["mu_XS7_r_XS8_WH"]
                        if P in ["qqH"]:
                            terms += ["mu_XS7_r_XS8_VBF"]
                    # Hack for ggH
                    if D in self.add_bbH and P == "ggH":
                        b2g = "CMS_R_bbH_ggH_%s_%dTeV[%g]" % (D, E, 0.01)
                        bbs = "CMS_bbH_scaler_%dTeV" % E
                        self.modelBuilder.factory_('expr::ggH_bbH_sum_%s_%dTeV("1+@0*@1",%s,%s)' % (D, E, b2g, bbs))
                        terms += ["ggH_bbH_sum_%s_%dTeV" % (D, E)]

                    self.modelBuilder.factory_("prod::scaling_%s_%s_%dTeV(%s)" % (P, D, E, ",".join(terms)))
                    self.modelBuilder.out.function("scaling_%s_%s_%dTeV" % (P, D, E)).Print("")

    def getHiggsSignalYieldScale(self, production, decay, energy):
        return f"scaling_{production}_{decay}_{energy}"


class Kappas(LHCHCGBaseModel):
    "assume the SM coupling but leave the Higgs mass to float"

    def __init__(
        self,
        resolved=True,
        BRU=True,
        addInvisible=False,
        addUndet=False,
        addWidth=False,
        addKappaC=False,
        custodial=False,
    ):
        LHCHCGBaseModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it
        self.doBRU = BRU
        self.resolved = resolved
        self.addInvisible = addInvisible
        self.addUndet = addUndet
        self.addKappaC = addKappaC
        self.addWidth = addWidth
        self.custodial = custodial

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("BRU="):
                self.doBRU = po.replace("BRU=", "") in [
                    "yes",
                    "1",
                    "Yes",
                    "True",
                    "true",
                ]
        print("BR uncertainties in partial widths: %s " % self.doBRU)

    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        if not self.custodial:
            self.modelBuilder.doVar("kappa_W[1,0.0,2.0]")
            self.modelBuilder.doVar("kappa_Z[1,0.0,2.0]")
            self.kappa_W = "kappa_W"
            self.kappa_Z = "kappa_Z"
        else:
            self.modelBuilder.doVar("kappa_V[1,0.0,2.0]")
            self.kappa_W = "kappa_V"
            self.kappa_Z = "kappa_V"
        kappa_b = "kappa_b"
        if self.addWidth:
            kappa_b = "c7_Gscal_tot"
            self.modelBuilder.doVar("c7_Gscal_tot[1,0.5,2.0]")
        else:
            self.modelBuilder.doVar("kappa_b[1,0.0,3.0]")
        self.modelBuilder.doVar("kappa_tau[1,0.0,3.0]")
        self.modelBuilder.doVar("kappa_mu[1,0.0,5.0]")
        # self.modelBuilder.factory_("expr::kappa_mu_expr(\"@0*@1+(1-@0)*@2\", CMS_use_kmu[0], kappa_mu, kappa_tau)")
        self.modelBuilder.doVar("kappa_t[1,0.0,4.0]")
        if not self.resolved:
            self.modelBuilder.doVar("kappa_g[1,0.0,2.0]")
            self.modelBuilder.doVar("kappa_gam[1,0.0,2.5]")
        self.modelBuilder.doVar("BRinv[0,0,1]")
        self.modelBuilder.doVar("BRundet[0,0,1]")
        if not self.addInvisible:
            self.modelBuilder.out.var("BRinv").setConstant(True)
        if not self.addUndet:
            self.modelBuilder.out.var("BRundet").setConstant(True)
        pois = self.kappa_W + "," + self.kappa_Z + ",kappa_tau,kappa_t," + kappa_b
        if not self.resolved:
            pois += ",kappa_g,kappa_gam"
        if self.addInvisible:
            pois += ",BRinv"
        if self.addUndet:
            pois += ",BRundet"
        if self.addKappaC:
            self.modelBuilder.doVar("kappa_c[1,0.0,10.0]")
            pois += ",kappa_c"
        if self.promote_hzg and not self.resolved:
            self.modelBuilder.doVar("kappa_Zgam[1.0,0.0,5.0]")
            pois += ",kappa_Zgam"
        self.doMH()
        self.modelBuilder.doSet("POI", pois)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        # SM BR
        for d in SM_HIGG_DECAYS + ["hss"]:
            self.SMH.makeBR(d)
        # BR uncertainties
        if self.doBRU:
            self.SMH.makePartialWidthUncertainties()
        else:
            for d in SM_HIGG_DECAYS:
                self.modelBuilder.factory_("HiggsDecayWidth_UncertaintyScaling_%s[1.0]" % d)
        # get VBF, tHq, tHW, ggZH cross section
        self.SMH.makeScaling("qqH", CW=self.kappa_W, CZ=self.kappa_Z)
        self.SMH.makeScaling("tHq", CW=self.kappa_W, Ctop="kappa_t")
        self.SMH.makeScaling("tHW", CW=self.kappa_W, Ctop="kappa_t")
        # resolve loops
        if self.resolved:
            if self.addKappaC:
                self.SMH.makeScaling("ggH", Cb="kappa_b", Ctop="kappa_t", Cc="kappa_c")
            else:
                self.SMH.makeScaling("ggH", Cb="kappa_b", Ctop="kappa_t", Cc="kappa_t")
            self.SMH.makeScaling("hgluglu", Cb="kappa_b", Ctop="kappa_t")
            self.SMH.makeScaling("hgg", Cb="kappa_b", Ctop="kappa_t", CW=self.kappa_W, Ctau="kappa_tau")
            self.SMH.makeScaling("hzg", Cb="kappa_b", Ctop="kappa_t", CW=self.kappa_W, Ctau="kappa_tau")
        else:
            self.modelBuilder.factory_('expr::Scaling_hgluglu("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_hgg("@0*@0", kappa_gam)')
            if self.promote_hzg:
                self.modelBuilder.factory_('expr::Scaling_hzg("@0*@0", kappa_Zgam)')
            else:
                self.modelBuilder.factory_('expr::Scaling_hzg("@0*@0", kappa_gam)')
            self.modelBuilder.factory_('expr::Scaling_ggH_7TeV("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_ggH_8TeV("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_ggH_13TeV("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_ggH_14TeV("@0*@0", kappa_g)')

        ## partial witdhs, normalized to the SM one
        kappa_mu_expr = "kappa_mu" if self.promote_hmm else "kappa_tau"
        self.modelBuilder.factory_('expr::c7_Gscal_Z("@0*@0*@1*@2", ' + self.kappa_Z + ", SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz)")
        self.modelBuilder.factory_('expr::c7_Gscal_W("@0*@0*@1*@2", ' + self.kappa_W + ", SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww)")
        if self.addKappaC:
            self.modelBuilder.factory_('expr::c7_Gscal_top("@0*@0 * @1*@2", kappa_c, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)')
        else:
            self.modelBuilder.factory_('expr::c7_Gscal_top("@0*@0 * @1*@2", kappa_t, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)')
        self.modelBuilder.factory_('expr::c7_Gscal_gluon("  @0  * @1 * @2", Scaling_hgluglu, SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu)')
        self.modelBuilder.factory_(
            'expr::c7_Gscal_gamma("@0*@1*@4 + @2*@3*@5",  Scaling_hgg, SM_BR_hgg, Scaling_hzg, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hgg, HiggsDecayWidth_UncertaintyScaling_hzg)'
        )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_tau("@0*@0*@1*@4+@2*@2*@3*@5", kappa_tau, SM_BR_htt, %s, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)'
            % kappa_mu_expr
        )
        # fix to have all BRs add up to unity
        self.modelBuilder.factory_("sum::c7_SMBRs(%s)" % (",".join("SM_BR_" + X for X in "hzz hww htt hmm hcc hbb hss hgluglu hgg hzg".split())))
        self.modelBuilder.out.function("c7_SMBRs").Print("")

        if self.addWidth:
            self.modelBuilder.factory_(
                'expr::c7_Gscal_bottom("@1*@8*(1-@0-@9)-(@2+@3+@4+@5+@6+@7)", BRinv, c7_Gscal_tot, c7_Gscal_W, c7_Gscal_Z, c7_Gscal_top, c7_Gscal_tau, c7_Gscal_gluon, c7_Gscal_gamma, c7_SMBRs, BRundet)'
            )
            self.modelBuilder.factory_('expr::kappa_b("sqrt(@0/(@1*@3+@2))", c7_Gscal_bottom, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb)')

        else:
            self.modelBuilder.factory_('expr::c7_Gscal_bottom("@0*@0 * (@1*@3+@2)", kappa_b, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb)')

        ## total witdh, normalized to the SM one
        if not self.addWidth:
            self.modelBuilder.factory_(
                'expr::c7_Gscal_tot("(@1+@2+@3+@4+@5+@6+@7)/@8/(1-@0-@9)", BRinv, c7_Gscal_Z, c7_Gscal_W, c7_Gscal_tau, c7_Gscal_top, c7_Gscal_bottom, c7_Gscal_gluon, c7_Gscal_gamma, c7_SMBRs, BRundet)'
            )

        self.SMH.makeScaling("ggZH", CZ=self.kappa_Z, Ctop="kappa_t", Cb="kappa_b")

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM) / (total/total_SM)
        self.modelBuilder.factory_('expr::c7_BRscal_hww("@0*@0*@2/@1", ' + self.kappa_W + ", c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hww)")
        self.modelBuilder.factory_('expr::c7_BRscal_hzz("@0*@0*@2/@1", ' + self.kappa_Z + ", c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzz)")
        self.modelBuilder.factory_('expr::c7_BRscal_htt("@0*@0*@2/@1", kappa_tau, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_htt)')
        self.modelBuilder.factory_('expr::c7_BRscal_hmm("@0*@0*@2/@1", %s, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hmm)' % kappa_mu_expr)
        self.modelBuilder.factory_('expr::c7_BRscal_hbb("@0*@0*@2/@1", kappa_b, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hbb)')
        if self.addKappaC:
            self.modelBuilder.factory_('expr::c7_BRscal_hcc("@0*@0*@2/@1", kappa_c, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc)')
        else:
            self.modelBuilder.factory_('expr::c7_BRscal_hcc("@0*@0*@2/@1", kappa_t, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc)')
        self.modelBuilder.factory_('expr::c7_BRscal_hgg("@0*@2/@1", Scaling_hgg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgg)')
        self.modelBuilder.factory_('expr::c7_BRscal_hzg("@0*@2/@1", Scaling_hzg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzg)')
        self.modelBuilder.factory_('expr::c7_BRscal_hgluglu("@0*@2/@1", Scaling_hgluglu, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgluglu)')

        self.modelBuilder.factory_('expr::c7_BRscal_hinv("@0", BRinv)')

    def getHiggsSignalYieldScale(self, production, decay, energy):
        name = f"c7_XSBRscal_{production}_{decay}_{energy}"
        if self.modelBuilder.out.function(name) == None:
            if production in ["ggH", "qqH", "ggZH", "tHq", "tHW"]:
                XSscal = ("@0", f"Scaling_{production}_{energy}")
            elif production == "WH":
                XSscal = ("@0*@0", self.kappa_W)
            elif production == "ZH":
                XSscal = ("@0*@0", self.kappa_Z)
            elif production == "ttH":
                XSscal = ("@0*@0", "kappa_t")
            elif production == "bbH":
                XSscal = ("@0*@0", "kappa_b")
            else:
                raise RuntimeError("Production %s not supported" % production)
            BRscal = decay
            if not self.modelBuilder.out.function("c7_BRscal_" + BRscal):
                raise RuntimeError("Decay mode %s not supported" % decay)
            if decay == "hss":
                BRscal = "hbb"
            if production == "ggH" and (decay in self.add_bbH) and energy in ["7TeV", "8TeV", "13TeV", "14TeV"]:
                b2g = f"CMS_R_bbH_ggH_{decay}_{energy}[{0.01:g}]"
                b2gs = "CMS_bbH_scaler_%s" % energy
                self.modelBuilder.factory_(f'expr::{name}("({XSscal[0]} + @1*@1*@2*@3)*@4", {XSscal[1]}, kappa_b, {b2g}, {b2gs}, c7_BRscal_{BRscal})')
            else:
                self.modelBuilder.factory_(f'expr::{name}("{XSscal[0]}*@1", {XSscal[1]}, c7_BRscal_{BRscal})')
            print("[LHC-HCG Kappas]", name, production, decay, energy, ": ", end=" ")
            self.modelBuilder.out.function(name).Print("")
        return name


class Lambdas(LHCHCGBaseModel):
    def __init__(self, BRU=True):
        LHCHCGBaseModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it
        self.doBRU = BRU

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("BRU="):
                self.doBRU = po.replace("BRU=", "") in [
                    "yes",
                    "1",
                    "Yes",
                    "True",
                    "true",
                ]
        print("BR uncertainties in partial widths: %s " % self.doBRU)

    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.doMH()
        self.modelBuilder.doVar("lambda_WZ[1,0.0,2.0]")
        self.modelBuilder.doVar("lambda_Zg[1,0.0,4.0]")
        self.modelBuilder.doVar("lambda_bZ[1,0.0,4.0]")
        self.modelBuilder.doVar("lambda_gamZ[1,0.0,2.0]")
        self.modelBuilder.doVar("lambda_tauZ[1,0.0,4.0]")
        self.modelBuilder.doVar("lambda_muZ[1,0.0,4.0]")
        self.modelBuilder.doVar("lambda_tg[1,0.0,4.0]")
        self.modelBuilder.doVar("kappa_gZ[1,0.0,3.0]")
        self.modelBuilder.doSet(
            "POI",
            "lambda_WZ,lambda_Zg,lambda_bZ,lambda_gamZ,lambda_tauZ,lambda_muZ,lambda_tg,kappa_gZ",
        )
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]), float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar(f"MH[{self.mHRange[0]},{self.mHRange[1]}]")
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        # BR uncertainties
        if self.doBRU:
            self.SMH.makePartialWidthUncertainties()

        self.modelBuilder.doVar("BRinv[0,0,1]")
        self.modelBuilder.out.var("BRinv").setConstant()

        self.modelBuilder.doVar("lambda_one[1]")

        self.modelBuilder.factory_('expr::lambda_tZ("@0/@1",lambda_tg,lambda_Zg)')
        self.SMH.makeScaling("qqH", CW="lambda_WZ", CZ="lambda_one")
        self.SMH.makeScaling("tHq", CW="lambda_WZ", Ctop="lambda_tZ")
        self.SMH.makeScaling("tHW", CW="lambda_WZ", Ctop="lambda_tZ")
        self.SMH.makeScaling("ggZH", CZ="lambda_one", Ctop="lambda_tZ", Cb="lambda_bZ")
        self.SMH.makeScaling("hzg", Cb="lambda_bZ", Ctop="lambda_tZ", CW="lambda_WZ", Ctau="lambda_tauZ")
        self.modelBuilder.factory_('expr::lambda_gZ("1/@0",lambda_Zg)')
        self.modelBuilder.factory_('expr::sqrt_zgamma("sqrt(@0)",Scaling_hzg)')

        for E in "7TeV", "8TeV", "13TeV", "14TeV":
            for P in "qqH", "tHq", "tHW", "ggZH":
                self.modelBuilder.factory_(f'expr::PW_XSscal_{P}_{E}("@0*@1*@1",Scaling_{P}_{E},lambda_Zg)')
            self.modelBuilder.factory_('expr::PW_XSscal_WH_%s("@0*@0*@1*@1",lambda_Zg,lambda_WZ)' % E)
            self.modelBuilder.factory_('expr::PW_XSscal_ZH_%s("@0*@0",lambda_Zg)' % E)
            self.modelBuilder.factory_('expr::PW_XSscal_ttH_%s("@0*@0",lambda_tg)' % E)
            self.modelBuilder.factory_('expr::PW_XSscal_bbH_%s("@0*@0*@1*@1",lambda_Zg,lambda_bZ)' % E)
        self.decayMap_ = {
            "hinv": "BRinv",
            "hww": "lambda_WZ",
            "hzz": "lambda_one",
            "hgg": "lambda_gamZ",
            "hbb": "lambda_bZ",
            "htt": "lambda_tauZ",
            "hmm": "lambda_muZ",
            "hcc": "lambda_tZ",  # charm scales as top
            "hgluglu": "lambda_gZ",  # glu scales as 1/Zgky
            "hzg": "lambda_gamZ",  # fancier option: 'sqrt_zgamma',
            #'hss' : 'lambda_bZ', # strange scales as bottom # not used
        }
        if not self.promote_hmm:
            self.decayMap_["hmm"] = "lambda_tauZ"

    def getHiggsSignalYieldScale(self, production, decay, energy):
        name = f"c7_XSBRscal_{production}_{decay}_{energy}"
        if self.modelBuilder.out.function(name):
            return name
        dscale = self.decayMap_[decay]
        if self.doBRU:
            name += "_noBRU"
        if production == "ggH":
            if decay in self.add_bbH:
                b2g = f"CMS_R_bbH_ggH_{decay}_{energy}[{0.01:g}]"
                b2gs = "CMS_bbH_scaler_%s" % energy
                self.modelBuilder.factory_(f'expr::{name}("@0*@0*@1*@1*(1+@2*@3*@4*@4*@5*@5)",kappa_gZ,{dscale},{b2g},{b2gs},lambda_bZ,lambda_Zg)')
            else:
                self.modelBuilder.factory_(f'expr::{name}("@0*@0*   @1*@1",kappa_gZ,{dscale})')
        else:
            self.modelBuilder.factory_(f'expr::{name}("@0*@0*@1*@2*@2",kappa_gZ,PW_XSscal_{production}_{energy},{dscale})')
        if self.doBRU:
            name = name.replace("_noBRU", "")
            if decay == "hzz":
                self.modelBuilder.factory_("prod::{}({}_noBRU, HiggsDecayWidth_UncertaintyScaling_{})".format(name, name, "hzz"))
            else:
                self.modelBuilder.factory_(
                    'expr::%s("@0*(@1/@2)", %s_noBRU, HiggsDecayWidth_UncertaintyScaling_%s, HiggsDecayWidth_UncertaintyScaling_%s)'
                    % (name, name, decay, "hzz")
                )
        print("[LHC-HCG Lambdas]", name, production, decay, energy, ": ", end=" ")
        self.modelBuilder.out.function(name).Print("")
        return name


class KappaVKappaF(LHCHCGBaseModel):
    "assume the SM coupling but let the Higgs mass to float"

    def __init__(self, BRU=True, floatbrinv=False):
        LHCHCGBaseModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it
        self.doBRU = BRU
        self.floatbrinv = floatbrinv

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("BRU="):
                self.doBRU = po.replace("BRU=", "") in [
                    "yes",
                    "1",
                    "Yes",
                    "True",
                    "true",
                ]
        print("BR uncertainties in partial widths: %s " % self.doBRU)

    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.modelBuilder.doVar("kappa_V[1,0.0,2.0]")
        self.modelBuilder.doVar("kappa_F[1,-2.0,2.0]")
        self.modelBuilder.doVar("BRinv[0,0,1]")
        self.modelBuilder.out.var("BRinv").setConstant()
        for d in ["WW", "ZZ", "gamgam", "bb", "tautau", "mumu", "inv"]:
            self.modelBuilder.doVar("kappa_V_%s[1,0.0,2.0]" % d)
            self.modelBuilder.doVar("kappa_F_%s[1,-2.0,2.0]" % d)
            self.modelBuilder.out.var("kappa_V_" + d).setConstant()
            self.modelBuilder.out.var("kappa_F_" + d).setConstant()
        pois = "kappa_V,kappa_F"
        if self.floatbrinv:
            self.modelBuilder.out.var("BRinv").setConstant(False)
            pois += ",BRinv"
        self.doMH()
        self.modelBuilder.doSet("POI", pois)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        # SM BR
        for d in SM_HIGG_DECAYS + ["hss"]:
            self.SMH.makeBR(d)
        # BR uncertainties
        if self.doBRU:
            self.SMH.makePartialWidthUncertainties()
        else:
            for d in SM_HIGG_DECAYS:
                self.modelBuilder.factory_("HiggsDecayWidth_UncertaintyScaling_%s[1.0]" % d)

        # fix to have all BRs add up to unity
        self.modelBuilder.factory_("sum::c7_SMBRs(%s)" % (",".join("SM_BR_" + X for X in "hzz hww htt hmm hcc hbb hss hgluglu hgg hzg".split())))
        self.modelBuilder.out.function("c7_SMBRs").Print("")

        for ds in ["WW", "ZZ", "gamgam", "bb", "tautau", "mumu", "inv"]:
            self.modelBuilder.factory_(f'expr::kVkV_{ds}("@0*@1", kappa_V,kappa_V_{ds})')
            self.modelBuilder.factory_(f'expr::kFkF_{ds}("@0*@1", kappa_F,kappa_F_{ds})')

            # get tHq, tHW, ggZH cross section
            self.SMH.makeScaling("tHq_" + ds, CW="kVkV_" + ds, Ctop="kFkF_" + ds)
            self.SMH.makeScaling("tHW_" + ds, CW="kVkV_" + ds, Ctop="kFkF_" + ds)
            self.SMH.makeScaling("ggZH_" + ds, CZ="kVkV_" + ds, Ctop="kFkF_" + ds, Cb="kFkF_" + ds)
            # resolve hgg, hzg loops
            self.SMH.makeScaling(
                "hgg_" + ds,
                Cb="kFkF_" + ds,
                Ctop="kFkF_" + ds,
                CW="kVkV_" + ds,
                Ctau="kFkF_" + ds,
            )
            self.SMH.makeScaling(
                "hzg_" + ds,
                Cb="kFkF_" + ds,
                Ctop="kFkF_" + ds,
                CW="kVkV_" + ds,
                Ctau="kFkF_" + ds,
            )
            ## partial witdhs, normalized to the SM one
            self.modelBuilder.factory_(f'expr::c7_Gscal_Z_{ds}("@0*@0*@1*@2", kVkV_{ds}, SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz)')
            self.modelBuilder.factory_(f'expr::c7_Gscal_W_{ds}("@0*@0*@1*@2", kVkV_{ds}, SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww)')
            self.modelBuilder.factory_(
                'expr::c7_Gscal_tau_%s("@0*@0*@1*@3+@0*@0*@2*@4", kFkF_%s, SM_BR_htt, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)'
                % (ds, ds)
            )
            self.modelBuilder.factory_(f'expr::c7_Gscal_top_{ds}("@0*@0 * @1*@2", kFkF_{ds}, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)')
            self.modelBuilder.factory_(
                f'expr::c7_Gscal_bottom_{ds}("@0*@0 * (@1*@3+@2)", kFkF_{ds}, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb)'
            )
            self.modelBuilder.factory_(f'expr::c7_Gscal_gluon_{ds}("  @0*@0*@1*@2", kFkF_{ds}, SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu)')
            self.modelBuilder.factory_(
                'expr::c7_Gscal_gamma_%s("@0*@1*@4 + @2*@3*@5",  Scaling_hgg_%s, SM_BR_hgg, Scaling_hzg_%s, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hgg, HiggsDecayWidth_UncertaintyScaling_hzg)'
                % (ds, ds, ds)
            )

            ## total witdh, normalized to the SM one
            self.modelBuilder.factory_(
                'expr::c7_Gscal_tot_%s("(@0+@1+@2+@3+@4+@5+@6)/@7/(1-@8)", c7_Gscal_Z_%s, c7_Gscal_W_%s, c7_Gscal_tau_%s, c7_Gscal_top_%s, c7_Gscal_bottom_%s, c7_Gscal_gluon_%s, c7_Gscal_gamma_%s, c7_SMBRs, BRinv)'
                % (ds, ds, ds, ds, ds, ds, ds, ds)
            )

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM) / (total/total_SM)
        for d in ["hww", "hzz"]:
            self.modelBuilder.factory_(
                'expr::c7_BRscal_%s("@0*@0*@2/@1", kVkV_%s, c7_Gscal_tot_%s, HiggsDecayWidth_UncertaintyScaling_%s)'
                % (d, CMS_to_LHCHCG_DecSimple[d], CMS_to_LHCHCG_DecSimple[d], d)
            )
        for d in ["htt", "hmm", "hbb", "hcc", "hgluglu"]:
            self.modelBuilder.factory_(
                'expr::c7_BRscal_%s("@0*@0*@2/@1", kFkF_%s, c7_Gscal_tot_%s, HiggsDecayWidth_UncertaintyScaling_%s)'
                % (d, CMS_to_LHCHCG_DecSimple[d], CMS_to_LHCHCG_DecSimple[d], d)
            )
        for d in ["hgg", "hzg"]:
            self.modelBuilder.factory_(
                'expr::c7_BRscal_%s("@0*@2/@1", Scaling_%s_%s, c7_Gscal_tot_%s, HiggsDecayWidth_UncertaintyScaling_%s)'
                % (d, d, CMS_to_LHCHCG_DecSimple[d], CMS_to_LHCHCG_DecSimple[d], d)
            )

        # H->invisible scaling
        self.modelBuilder.factory_('expr::c7_BRscal_hinv("@0", BRinv)')

    def getHiggsSignalYieldScale(self, production, decay, energy):
        name = f"c7_XSBRscal_{production}_{decay}_{energy}"
        if self.modelBuilder.out.function(name) == None:
            if production in ["ggZH", "tHq", "tHW"]:
                XSscal = (
                    "@0",
                    f"Scaling_{production}_{CMS_to_LHCHCG_DecSimple[decay]}_{energy}",
                )
            elif production in ["ggH", "ttH", "bbH"]:
                XSscal = ("@0*@0", "kFkF_" + CMS_to_LHCHCG_DecSimple[decay])
            elif production in ["qqH", "WH", "ZH"]:
                XSscal = ("@0*@0", "kVkV_" + CMS_to_LHCHCG_DecSimple[decay])
            else:
                raise RuntimeError("Production %s not supported" % production)
            BRscal = decay
            if decay == "hss":
                BRscal = "hbb"
            if not self.modelBuilder.out.function("c7_BRscal_" + BRscal):
                raise RuntimeError("Decay mode %s not supported" % decay)
            if production == "ggH" and (decay in self.add_bbH) and energy in ["7TeV", "8TeV", "13TeV", "14TeV"]:
                b2g = f"CMS_R_bbH_ggH_{decay}_{energy}[{0.01:g}]"
                b2gs = "CMS_bbH_scaler_%s" % energy
                self.modelBuilder.factory_(
                    'expr::%s("(%s + @1*@1*@2*@3)*@4", %s, kFkF_%s, %s, %s, c7_BRscal_%s)'
                    % (
                        name,
                        XSscal[0],
                        XSscal[1],
                        CMS_to_LHCHCG_DecSimple[decay],
                        b2g,
                        b2gs,
                        BRscal,
                    )
                )
            else:
                self.modelBuilder.factory_(f'expr::{name}("{XSscal[0]}*@1", {XSscal[1]}, c7_BRscal_{BRscal})')
            print("[LHC-HCG Kappas]", name, production, decay, energy, ": ", end=" ")
            self.modelBuilder.out.function(name).Print("")
        return name


class LambdasReduced(LHCHCGBaseModel):
    def __init__(self, BRU=True, model="", addInvisible=False, flipped=False):
        LHCHCGBaseModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it
        self.doBRU = BRU
        self.model = model
        self.addInvisible = addInvisible
        self.flipped = flipped

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("BRU="):
                self.doBRU = po.replace("BRU=", "") in [
                    "yes",
                    "1",
                    "Yes",
                    "True",
                    "true",
                ]
        print("BR uncertainties in partial widths: %s " % self.doBRU)

    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.doMH()
        self.modelBuilder.doVar("lambda_FV[1,0.0,2.0]")
        self.modelBuilder.doVar("kappa_VV[1,0.0,2.0]")
        self.modelBuilder.doVar("lambda_du[1,0.0,2.0]")
        self.modelBuilder.doVar("lambda_Vu[1,0.0,2.0]")
        self.modelBuilder.doVar("kappa_uu[1,0.0,2.0]")
        self.modelBuilder.doVar("lambda_lq[1,0.0,2.0]")
        self.modelBuilder.doVar("lambda_Vq[1,0.0,2.0]")
        self.modelBuilder.doVar("kappa_qq[1,0.0,2.0]")

        self.modelBuilder.doVar("BRinv[0,0,1]")
        if not self.addInvisible:
            self.modelBuilder.out.var("BRinv").setConstant(True)

        POIset = "lambda_FV,kappa_VV,lambda_du,lambda_Vu,kappa_uu,lambda_lq,lambda_Vq,kappa_qq"
        if self.model == "ldu":
            self.modelBuilder.out.var("lambda_du").setConstant(False)
            self.modelBuilder.out.var("kappa_uu").setConstant(False)
            self.modelBuilder.out.var("lambda_Vu").setConstant(False)

            self.modelBuilder.out.var("lambda_FV").setConstant(True)
            self.modelBuilder.out.var("kappa_VV").setConstant(True)
            self.modelBuilder.out.var("lambda_lq").setConstant(True)
            self.modelBuilder.out.var("lambda_Vq").setConstant(True)
            self.modelBuilder.out.var("kappa_qq").setConstant(True)

            POIset = "lambda_du,lambda_Vu,kappa_uu"

        elif self.model == "llq":
            self.modelBuilder.out.var("lambda_lq").setConstant(False)
            self.modelBuilder.out.var("lambda_Vq").setConstant(False)
            self.modelBuilder.out.var("kappa_qq").setConstant(False)

            self.modelBuilder.out.var("lambda_du").setConstant(True)
            self.modelBuilder.out.var("kappa_uu").setConstant(True)
            self.modelBuilder.out.var("lambda_Vu").setConstant(True)
            self.modelBuilder.out.var("lambda_FV").setConstant(True)
            self.modelBuilder.out.var("kappa_VV").setConstant(True)

            POIset = "lambda_lq,lambda_Vq,kappa_qq"

        elif self.model == "lfv":
            self.modelBuilder.out.var("lambda_FV").setConstant(False)
            self.modelBuilder.out.var("kappa_VV").setConstant(False)

            self.modelBuilder.out.var("lambda_lq").setConstant(True)
            self.modelBuilder.out.var("lambda_Vq").setConstant(True)
            self.modelBuilder.out.var("kappa_qq").setConstant(True)
            self.modelBuilder.out.var("lambda_du").setConstant(True)
            self.modelBuilder.out.var("kappa_uu").setConstant(True)
            self.modelBuilder.out.var("lambda_Vu").setConstant(True)

            POIset = "lambda_FV,kappa_VV"

        self.modelBuilder.doSet("POI", POIset)
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]), float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar(f"MH[{self.mHRange[0]},{self.mHRange[1]}]")
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        # BR uncertainties
        if self.doBRU:
            self.SMH.makePartialWidthUncertainties()

        if self.flipped:
            self.modelBuilder.factory_('expr::C_b("@0",lambda_FV)')
            self.modelBuilder.factory_('expr::C_top("@0*@1",lambda_du,lambda_FV)')
            self.modelBuilder.factory_('expr::C_tau("@0*@1*@2",lambda_lq,lambda_du,lambda_FV)')
            self.modelBuilder.factory_('expr::C_V("@0*@1",lambda_Vq,lambda_Vu)')
        else:
            self.modelBuilder.factory_('expr::C_b("@0*@1",lambda_du,lambda_FV)')
            self.modelBuilder.factory_('expr::C_top("@0",lambda_FV)')
            self.modelBuilder.factory_('expr::C_tau("@0*@1*@2",lambda_lq,lambda_du, lambda_FV)')
            self.modelBuilder.factory_('expr::C_V("@0*@1",lambda_Vq,lambda_Vu)')

        self.SMH.makeScaling("ggH", Cb="C_b", Ctop="C_top", Cc="C_top")
        self.SMH.makeScaling("tHq", CW="C_V", Ctop="C_top")
        self.SMH.makeScaling("tHW", CW="C_V", Ctop="C_top")
        self.SMH.makeScaling("ggZH", CZ="C_V", Ctop="C_top", Cb="C_b")
        self.SMH.makeScaling("hgg", Cb="C_b", Ctop="C_top", CW="C_V", Ctau="C_tau")
        self.SMH.makeScaling("hzg", Cb="C_b", Ctop="C_top", CW="C_V", Ctau="C_tau")
        self.SMH.makeScaling("hgluglu", Cb="C_b", Ctop="C_top")

        # H->invisible scaling
        self.modelBuilder.factory_('expr::Scaling_hinv("@0", BRinv)')

        for E in "7TeV", "8TeV", "13TeV", "14TeV":
            for P in "ggH", "tHq", "tHW", "ggZH":
                self.modelBuilder.factory_(f'expr::PW_XSscal_{P}_{E}("@0*@1*@1*@2*@2*@3*@3",Scaling_{P}_{E},kappa_qq, kappa_uu, kappa_VV)')
            for P in "qqH", "WH", "ZH":
                self.modelBuilder.factory_(f'expr::PW_XSscal_{P}_{E}("@0*@0*@1*@1*@2*@2*@3*@3*@4*@4", kappa_qq, lambda_Vq, kappa_uu, lambda_Vu, kappa_VV)')
            if self.flipped:
                self.modelBuilder.factory_('expr::PW_XSscal_bbH_%s("@0*@0*@1*@1*@2*@2*@3*@3",kappa_qq,kappa_uu,kappa_VV,lambda_FV)' % E)
                self.modelBuilder.factory_('expr::PW_XSscal_ttH_%s("@0*@0*@1*@1*@2*@2*@3*@3*@4*@4",kappa_qq,kappa_uu,lambda_du,kappa_VV,lambda_FV)' % E)
            else:
                self.modelBuilder.factory_('expr::PW_XSscal_ttH_%s("@0*@0*@1*@1*@2*@2*@3*@3",kappa_qq,kappa_uu,kappa_VV,lambda_FV)' % E)
                self.modelBuilder.factory_('expr::PW_XSscal_bbH_%s("@0*@0*@1*@1*@2*@2*@3*@3*@4*@4",kappa_qq,kappa_uu,lambda_du,kappa_VV,lambda_FV)' % E)

        self.decayMap_ = {
            "hww": "C_V",
            "hzz": "C_V",
            "hgg": "Scaling_hgg",
            "hbb": "C_b",
            "htt": "C_tau",
            "hmm": "C_tau",
            "hcc": "C_top",  # charm scales as top
            "hgluglu": "Scaling_hgluglu",
            "hzg": "Scaling_hzg",
            "hinv": "Scaling_hinv",
            #'hss' : 'C_b', # strange scales as bottom # not used
        }

    def getHiggsSignalYieldScale(self, production, decay, energy):
        name = f"c7_XSBRscal_{production}_{decay}_{energy}"
        if self.modelBuilder.out.function(name):
            return name
        dscale = self.decayMap_[decay]
        if self.doBRU:
            name += "_noBRU"
        if production == "ggH":
            if decay in self.add_bbH:
                b2g = f"CMS_R_bbH_ggH_{decay}_{energy}[{0.01:g}]"
                b2gs = "CMS_bbH_scaler_%s" % energy
                if decay in ["hgg", "hgluglu", "hzg"]:
                    self.modelBuilder.factory_(f'expr::{name}("@0*@1*(1+@2*@3*@4)",PW_XSscal_ggH_{energy},{dscale},{b2g},{b2gs},PW_XSscal_bbH_{energy})')
                else:
                    self.modelBuilder.factory_(f'expr::{name}("@0*@1*@1*(1+@2*@3*@4)",PW_XSscal_ggH_{energy},{dscale},{b2g},{b2gs},PW_XSscal_bbH_{energy})')
            else:
                if decay in ["hgg", "hgluglu", "hzg"]:
                    self.modelBuilder.factory_(f'expr::{name}("@0*@1",PW_XSscal_ggH_{energy},{dscale})')
                else:
                    self.modelBuilder.factory_(f'expr::{name}("@0*@1*@1",PW_XSscal_ggH_{energy},{dscale})')
        else:
            if decay in ["hgg", "hgluglu", "hzg"]:
                self.modelBuilder.factory_(f'expr::{name}("@0*@1",PW_XSscal_{production}_{energy},{dscale})')
            else:
                self.modelBuilder.factory_(f'expr::{name}("@0*@1*@1",PW_XSscal_{production}_{energy},{dscale})')
        if self.doBRU:
            name = name.replace("_noBRU", "")
            if decay == "hzz":
                self.modelBuilder.factory_("prod::{}({}_noBRU, HiggsDecayWidth_UncertaintyScaling_{})".format(name, name, "hzz"))
            else:
                self.modelBuilder.factory_(
                    'expr::%s("@0*(@1/@2)", %s_noBRU, HiggsDecayWidth_UncertaintyScaling_%s, HiggsDecayWidth_UncertaintyScaling_%s)'
                    % (name, name, decay, "hzz")
                )
        print("[LHC-HCG Lambdas]", name, production, decay, energy, ": ", end=" ")
        self.modelBuilder.out.function(name).Print("")
        return name


class XSBRratiosAlternative(LHCHCGBaseModel):
    "Model with the ratios of cross sections and branching ratios"

    def __init__(self):
        LHCHCGBaseModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("poi="):
                self.POIs = po.replace("poi=", "")

    def doVar(self, x, constant=True):
        self.modelBuilder.doVar(x)
        vname = re.sub(r"\[.*", "", x)
        self.modelBuilder.out.var(vname).setConstant(constant)
        print(f"XSBRratios:: declaring {vname} as {x}, and set to constant")

    def doParametersOfInterest(self):
        """Create POI out of signal strength ratios and MH"""
        self.doMH()
        self.doVar("mu_XS7_r_XS8_ggF[1,0,5]")
        self.doVar("mu_XS7_r_XS8_VBF[1,0,5]")
        self.doVar("mu_XS7_r_XS8_WH[1,0,5]")
        self.doVar("mu_XS7_r_XS8_ZH[1,0,5]")
        self.doVar("mu_XS7_r_XS8_ttH[1,0,5]")
        self.modelBuilder.doVar("mu_XS_ggF_x_BR_WW[1,0,5]")
        self.modelBuilder.doVar("mu_XS_VBF_x_BR_tautau[1,0,5]")
        self.modelBuilder.doVar("mu_XS_WH_r_XS_VBF[1,0,5]")
        self.modelBuilder.doVar("mu_XS_ZH_r_XS_WH[1,0,5]")
        self.modelBuilder.doVar("mu_XS_ttH_r_XS_ggF[1,0,5]")
        self.POIs = "mu_XS_ggF_x_BR_WW,mu_XS_VBF_x_BR_tautau,mu_XS_ttH_r_XS_ggF,mu_XS_WH_r_XS_VBF,mu_XS_ZH_r_XS_WH"
        for X in ["ZZ", "tautau", "gamgam"]:
            self.modelBuilder.doVar("mu_BR_%s_r_BR_WW[1,0,5]" % X)
            self.POIs += ",mu_BR_%s_r_BR_WW" % (X)
        self.modelBuilder.doVar("mu_BR_bb_r_BR_tautau[1,0,5]")
        self.POIs += ",mu_BR_bb_r_BR_tautau"
        print("Default parameters of interest: ", self.POIs)
        self.modelBuilder.doSet("POI", self.POIs)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        for P in ALL_HIGGS_PROD:
            if P == "VH":
                continue  # skip aggregates
            for D in SM_HIGG_DECAYS:
                for E in 7, 8, 13:
                    terms = []
                    if P in ["ggH", "bbH", "ttH", "tHq", "tHW"] and CMS_to_LHCHCG_DecSimple[D] in ["WW", "ZZ", "gamgam"]:
                        terms = ["mu_XS_ggF_x_BR_WW"]
                        if CMS_to_LHCHCG_DecSimple[D] in ["ZZ", "gamgam"]:
                            terms += ["mu_BR_%s_r_BR_WW" % CMS_to_LHCHCG_DecSimple[D]]
                        if P in ["ttH", "tHq", "tHW"]:
                            terms += ["mu_XS_ttH_r_XS_ggF"]
                    elif P in ["ggH", "bbH", "ttH", "tHq", "tHW"] and CMS_to_LHCHCG_DecSimple[D] in ["tautau", "bb"]:
                        terms = ["mu_XS_ggF_x_BR_WW", "mu_BR_tautau_r_BR_WW"]
                        if CMS_to_LHCHCG_DecSimple[D] == "bb":
                            terms += ["mu_BR_bb_r_BR_tautau"]
                        if P in ["ttH", "tHq", "tHW"]:
                            terms += ["mu_XS_ttH_r_XS_ggF"]
                    elif P in ["qqH", "WH", "ZH", "ggZH"] and CMS_to_LHCHCG_DecSimple[D] in ["WW", "ZZ", "gamgam"]:
                        terms = ["mu_XS_VBF_x_BR_tautau", "1/BR_tautau_r_BR_WW"]
                        if CMS_to_LHCHCG_DecSimple[D] in ["ZZ", "gamgam"]:
                            terms += ["mu_BR_%s_r_BR_WW" % CMS_to_LHCHCG_DecSimple[D]]
                        if P in ["WH", "ZH", "ggZH"]:
                            terms += ["mu_XS_WH_r_XS_VBF"]
                            if P in ["ZH", "ggZH"]:
                                terms += ["mu_XS_ZH_r_XS_WH"]
                    elif P in ["qqH", "WH", "ZH", "ggZH"] and CMS_to_LHCHCG_DecSimple[D] in ["bb", "tautau"]:
                        terms = ["mu_XS_VBF_x_BR_tautau"]
                        if CMS_to_LHCHCG_DecSimple[D] in ["bb"]:
                            terms += ["mu_BR_bb_r_BR_tautau"]
                        if P in ["WH", "ZH", "ggZH"]:
                            terms += ["mu_XS_WH_r_XS_VBF"]
                            if P in ["ZH", "ggZH"]:
                                terms += ["mu_XS_ZH_r_XS_WH"]
                    if E == 7:
                        if P in ["ggH", "bbH"]:
                            terms += ["mu_XS7_r_XS8_ggF"]
                        if P in ["ttH", "tHq", "tHW"]:
                            terms += ["mu_XS7_r_XS8_ttH"]
                        if P in ["ZH", "ggZH"]:
                            terms += ["mu_XS7_r_XS8_ZH"]
                        if P in ["WH"]:
                            terms += ["mu_XS7_r_XS8_WH"]
                        if P in ["qqH"]:
                            terms += ["mu_XS7_r_XS8_VBF"]
                    # Hack for ggH
                    if D in self.add_bbH and P == "ggH":
                        b2g = "CMS_R_bbH_ggH_%s_%dTeV[%g]" % (D, E, 0.01)
                        bbs = "CMS_bbH_scaler_%dTeV" % E
                        self.modelBuilder.factory_('expr::ggH_bbH_sum_%s_%dTeV("1+@0*@1",%s,%s)' % (D, E, b2g, bbs))
                        terms += ["ggH_bbH_sum_%s_%dTeV" % (D, E)]

                    self.modelBuilder.factory_("prod::scaling_%s_%s_%dTeV(%s)" % (P, D, E, ",".join(terms)))
                    self.modelBuilder.out.function("scaling_%s_%s_%dTeV" % (P, D, E)).Print("")

    def getHiggsSignalYieldScale(self, production, decay, energy):
        return f"scaling_{production}_{decay}_{energy}"


class CommonMatrixModel(LHCHCGBaseModel):
    """Model for testing for degenerate states"""

    def __init__(self):
        LHCHCGBaseModel.__init__(self)  # not using 'super(x,self).__init__' since I don't understand it
        # set the list of POIs as empty
        self.POIs = []
        # set the list of missing ggH measurements as empty
        self.fixDecays = []

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("fixDecays="):
                self.fixDecays = po.replace("fixDecays=", "").split(",")
                for decay in self.fixDecays:
                    self.fixDecays[self.fixDecays.index(decay)] = CMS_to_LHCHCG_DecSimple[decay]
            if po.startswith("poi="):
                self.POIs = po.replace("poi=", "")

    def doVar(self, x, constant=True):
        self.modelBuilder.doVar(x)
        vname = re.sub(r"\[.*", "", x)
        self.modelBuilder.out.var(vname).setConstant(constant)
        print(f"DegenerateMatrixModel:: declaring {vname} as {x}, and set to constant")

    def doParametersOfInterest(self):
        """Create POI out of l_j and l_j_i, and mu_i"""
        for X in ["VBF", "WH", "ZH", "ttH"]:
            self.doVar("l_%s[1,0,10]" % X)
            self.POIs.append("l_%s" % X)
            for Y in ["gamgam", "WW", "ZZ", "tautau", "bb"]:
                self.doVar(f"l_{X}_{Y}[1,0,10]")
                self.POIs.append(f"l_{X}_{Y}")
        for D in ["gamgam", "WW", "ZZ", "tautau", "bb"]:
            self.modelBuilder.doVar("mu_%s[1,0,5]" % D)
            if D in self.fixDecays:  # If there are any missing ggH measurements, create the mus here and a POI out of it.
                print(
                    "It seems that you set signal strength ggH->"
                    + D
                    + " as missing. Take this into account when running HybridNew and using --redefineSignalPOIs!"
                )
                self.POIs.append("mu_%s" % D)
        print("The possible parameters of interest: ", self.POIs)
        self.modelBuilder.doSet("POI", ",".join(self.POIs))
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        for P in ALL_HIGGS_PROD:
            if P == "VH":
                continue  # skip aggregates
            for D in SM_HIGG_DECAYS:
                for E in 7, 8, 13:
                    terms = []
                    # Hack for ggH
                    if D in self.add_bbH and P == "ggH":
                        b2g = "CMS_R_bbH_ggH_%s_%dTeV[%g]" % (
                            CMS_to_LHCHCG_DecSimple[D],
                            E,
                            0.01,
                        )
                        bbs = "CMS_bbH_scaler_%dTeV" % E
                        ## FIXME should include the here also logNormal for QCDscale_bbH
                        self.modelBuilder.factory_(
                            'expr::scaler_%s_%dTeV("(@0+@1*@2)*@3",%d,%s,%s,mu_%s)'
                            % (
                                CMS_to_LHCHCG_DecSimple[D],
                                E,
                                1.0,
                                b2g,
                                bbs,
                                CMS_to_LHCHCG_DecSimple[D],
                            )
                        )
                        terms += [
                            "scaler_%s_%dTeV" % (CMS_to_LHCHCG_DecSimple[D], E),
                            "mu_" + CMS_to_LHCHCG_DecSimple[D],
                        ]
                    else:
                        if P in ["ggH", "bbH"]:
                            terms += ["mu_" + CMS_to_LHCHCG_DecSimple[D]]
                    # Parametrizations for rest of the rows
                    if P == "qqH":
                        terms += [
                            "mu_" + CMS_to_LHCHCG_DecSimple[D],
                            "l_" + CMS_to_LHCHCG_Prod[P],
                            f"l_{CMS_to_LHCHCG_Prod[P]}_{CMS_to_LHCHCG_DecSimple[D]}",
                        ]
                    if P in ["ZH", "ggZH"]:
                        # Will the production mode be now ZH, qqZH or what? See naming conventions... qqHz missing completely in ALL_HIGGS_PROD
                        terms += [
                            "mu_" + CMS_to_LHCHCG_DecSimple[D],
                            "l_ZH",
                            "l_ZH_%s" % CMS_to_LHCHCG_DecSimple[D],
                        ]
                    if P == "WH":
                        terms += [
                            "mu_" + CMS_to_LHCHCG_DecSimple[D],
                            "l_" + CMS_to_LHCHCG_Prod[P],
                            f"l_{CMS_to_LHCHCG_Prod[P]}_{CMS_to_LHCHCG_DecSimple[D]}",
                        ]
                    if P in ["ttH", "tHq", "tHW"]:
                        terms += [
                            "mu_" + CMS_to_LHCHCG_DecSimple[D],
                            "l_ttH",
                            "l_ttH_%s" % CMS_to_LHCHCG_DecSimple[D],
                        ]
                    self.modelBuilder.factory_("prod::scaling_%s_%s_%dTeV(%s)" % (P, D, E, ",".join(terms)))
                    self.modelBuilder.out.function("scaling_%s_%s_%dTeV" % (P, D, E)).Print("")

    def getHiggsSignalYieldScale(self, production, decay, energy):
        return f"scaling_{production}_{decay}_{energy}"


A1 = SignalStrengths()
A2 = SignalStrengthRatios()
B1 = XSBRratios("WW")
B1ZZ = XSBRratios("ZZ")
B2 = XSBRratiosAlternative()
K1 = Kappas(resolved=True)
K2 = Kappas(resolved=False)
K2Width = Kappas(resolved=False, addWidth=True)
K2Inv = Kappas(resolved=False, addInvisible=True, addUndet=False)
K2InvC = Kappas(resolved=False, addInvisible=True, addUndet=False, addWidth=False, addKappaC=True)
K2InvWidth = Kappas(resolved=False, addInvisible=True, addUndet=False, addWidth=True)
K2Undet = Kappas(resolved=False, addInvisible=True, addUndet=True)
K2UndetWidth = Kappas(resolved=False, addInvisible=True, addUndet=True, addWidth=True)
K3 = KappaVKappaF(floatbrinv=False)
K3Inv = KappaVKappaF(floatbrinv=True)
L1 = Lambdas()
L2 = LambdasReduced()
L2flipped = LambdasReduced(flipped=True)
D1 = CommonMatrixModel()
