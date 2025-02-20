/* includes */
#include <iostream>

#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomp.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_Fitters.h"
#include "DLM_HistoAnalysis.h"
#include "DLM_RootWrapper.h"
#include "CommonAnaFunctions.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TString.h"
#include "TObjString.h"

#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/Fitter.h"
#include "HFitInterface.h"

#include "cf_fitter.h"

/* defines  */
#define CUSP_WEIGHT 0.33
#define LAM_PL_FLAT 0.4969
#define CUTOFF_Q    200

#define POT_PP	    "AV18"
#define POT_PL	    "Chiral_Coupled_SPD"
#define POT_PS	    "DG_NLO19"
#define SOURCE	    "Gauss"
#define SOURCE_RSM  "McLevy_ResoTM"

#define out(X) std::cout << X << std::endl

double bslfun(double *mom, double *par)
{
    double &k = *mom;

    double &N = par[0];
    double &a = par[1];
    double &b = par[2];
    double &c = par[3];
    double &d = par[4];

    return N*(1. + a*k + b*pow(k, 2) + c*pow(k, 3.) + d*pow(k, 4.));
}

class plFitFunctionSimple
{
    public:
	DLM_CkDecomposition *mydecomp;
	plFitFunctionSimple(DLM_CkDecomposition *comp) : mydecomp (comp) {}

	double operator() (double *mom, double *par)
	{
	    double &k = *mom;
	    double &N = par[0];
	    double &a = par[1];
	    double &b = par[2];
	    double &c = par[3];
	    double &d = par[4];
	    double &radius = par[5];

	    mydecomp->GetCk()->SetSourcePar(0, radius);
	    mydecomp->Update(false, false);

	    return N*(1. + a*k + b*k*k + c*k*k*k + d*k*k*k*k) * (mydecomp->EvalCk(k));
	}
};

class FitFunSimpleAvg
{
    public:
	DLM_CkDecomposition *mydecomp;
	TH1F *cf;
	TH1F *me;

	FitFunSimpleAvg(TH1F *hData, TH1F *hME, DLM_CkDecomposition *comp)
	{
	    mydecomp = comp;
	    cf = static_cast<TH1F*>(hData->Clone("correlation_function"));
	    me = static_cast<TH1F*>(hME->Clone("mixed_event"));
	}

	FitFunSimpleAvg(TH1F **hData, TH1F **hMe, DLM_CkDecomposition *comp)
	{
	    mydecomp = comp;
	    cf = static_cast<TH1F*>((*hData)->Clone("correlation_function"));
	    me = static_cast<TH1F*>((*hMe)->Clone("mixed_event"));
	}

	double operator() (double *mom, double *par)
	{
	    double &k = *mom;
	    double &N = par[0];
	    double &a = par[1];
	    double &b = par[2];
	    double &c = par[3];
	    double &d = par[4];
	    double &radius = par[5];

	    mydecomp->GetCk()->SetSourcePar(0, radius);
	    mydecomp->Update(true, true);

	    double cf_avg = 0;
	    double me_sum = 0;
	    int bin = cf->GetXaxis()->FindBin(k);

	    double value_low = cf->GetXaxis()->GetBinLowEdge(bin);
	    double value_up  = cf->GetXaxis()->GetBinLowEdge(bin + 1);
	    int me_bin_low = me->GetXaxis()->FindBin(value_low);
	    int me_bin_up  = me->GetXaxis()->FindBin(value_up);

	    for (int nbin = me_bin_low; nbin < me_bin_up; ++nbin)
	    {
		cf_avg += me->GetBinContent(nbin) * mydecomp->EvalCk(me->GetBinCenter(nbin));
		me_sum += me->GetBinContent(nbin);
	    }
	    double ck_eval = (me_sum < 1e-6)? 0 : cf_avg / me_sum;

	    return N*(1. + a*k + b*pow(k, 2) + c*pow(k, 3.) + d*pow(k, 4.)) * ck_eval;
	}
};

class GenuineFun
{
    public:
	DLM_CkDecomposition *mydecomp;

	GenuineFun(DLM_CkDecomposition *comp) { mydecomp = comp; }
	void set_decomp(DLM_CkDecomposition *comp) { mydecomp = comp; }

	double operator() (double *mom, double *par)
	{
	    double &k = *mom;
	    double &N = par[0];
	    return N*(mydecomp->GetCk()->Eval(k));
	}
};

class GenuineFunSmear
{
    public:
	DLM_CkDecomposition *mydecomp;

	GenuineFunSmear(DLM_CkDecomposition *comp) { mydecomp = comp; }
	void set_decomp(DLM_CkDecomposition *comp) { mydecomp = comp; }

	double operator() (double *mom, double *par)
	{
	    double &k = *mom;
	    double &lambda = par[0];

	    return lambda * (mydecomp->EvalSmearedMain(k)) + (1. - lambda);
	}
};

class GlobalChi2
{
    public:
	int *fit_pars_pp;
	int *fit_pars_aa;

	ROOT::Math::IMultiGenFunction *fChi2_1;
	ROOT::Math::IMultiGenFunction *fChi2_2;

	GlobalChi2(ROOT::Math::IMultiGenFunction &f1, ROOT::Math::IMultiGenFunction &f2, int *par_pp, int *par_aa)
	{
	    fit_pars_pp = par_pp;
	    fit_pars_aa = par_aa;

	    fChi2_1 = &f1;
	    fChi2_2 = &f2;
	}

	double operator() (const double *par) const
	{
	    double par1[6], par2[6];
	    for (int i = 0; i < 6; ++i)
	    {
		par1[i] = par[fit_pars_pp[i]];
		par2[i] = par[fit_pars_aa[i]];
	    }

	    return (*fChi2_1)(par1) + (*fChi2_2)(par2);
	}
};

void print_info_2(VAR *var, VAR_FMR *fmr, VAR_FR *fr)
{
    if (var->sample)
    {
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │ Random Sample │\e[0m   %i\n\e[0m", var->sample);
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");
    }

    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  System       │\e[0m   %i  (%s)\n\e[0m", var->system, SYSTEM_STR[var->system].Data());
    printf("\e[1;34m  │  Charge       │\e[0m   %i  (%s)\n\e[0m", var->charge, CHARGE_STR[var->charge].Data());
    printf("\e[1;34m  │  Bin mt       │\e[0m   %i\n\e[0m", var->mt);
    printf("\e[1;34m  │  Bin mult     │\e[0m   %i\n\e[0m", var->mult);
    printf("\e[1;34m  ├───────────────┤\e[0m\n");
    printf("\e[1;34m  │  LambdaPar    │\e[0m   %i  (%s)\n\e[0m", var->lam, (var->lam == 1)? "-10%" : (var->lam == 2)? "+10%" : "default");
    printf("\e[1;34m  │  Femtorange   │\e[0m   %i  (%.0f, %.0f)\n\e[0m", var->fmr, fmr->Min, fmr->Max);
    printf("\e[1;34m  │  Fitrange     │\e[0m   %i  (%.0f, %.0f)\n\e[0m", var->fr, fr->Min, fr->Max);
    printf("\e[1;34m  │  Baseline     │\e[0m   %i  (%s)\n\e[0m", var->bsl, (var->bsl)? "Pol2" : "Pol3");
    printf("\e[1;34m  │  Smearing     │\e[0m   %i  (%s)\n\e[0m", var->smear, (var->smear == 1)? "Run3" : (var->smear == 2)? "Run2" : "None");

    if (var->rsm)
    {
    printf("\e[1;34m  ├───────────────┤\e[0m\n");
    printf("\e[1;34m  │  RSM Frac     │\e[0m   %i  (%.4f)\n\e[0m", var->frac, RSM_FRAC[var->system][var->frac]);
    printf("\e[1;34m  │  RSM Mass     │\e[0m   %i  (%.1f)\n\e[0m", var->mass, RSM_MASS[var->system][var->mass]);
    }
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");
}

void get_klambda_histos(TString filename, VAR *var, DLM_Histo<double> *lambda_pars[])
{
    // k* dependent λ parameters
    TString klambda_filename = filename;
    auto file = std::make_unique<TFile>(klambda_filename, "read");
    auto lambdas = std::make_unique<TDirectory*[]>(5);
    lambdas[0] = (TDirectory*) file->GetDirectory("genuine");
    lambdas[1] = (TDirectory*) file->GetDirectory("FD_lambda");
    lambdas[2] = (TDirectory*) file->GetDirectory("FD_sigma");
    lambdas[3] = (TDirectory*) file->GetDirectory("flat");
    lambdas[4] = (TDirectory*) file->GetDirectory("fake");

    double scale = 0.1;
    for (size_t nlam = 0; nlam < 5; ++nlam)
    {
	TString name = Form("Mult_%i/Mt_%i/Lambda", var->mult, var->mt - 1);
	auto th1_klambda = (TH1F*) lambdas[nlam]->Get(name);
	th1_klambda->SetDirectory(0);

	if	(var->lam == 1) th1_klambda->Scale(1. - scale);
	else if (var->lam == 2) th1_klambda->Scale(1. + scale);

	lambda_pars[nlam] = Convert_TH1F_DoubleDlmHisto(th1_klambda);
    }
    file->Close();
}

void get_klambda_histos(TString filename, VAR *var, TH1F *lambda_pars[])
{
    // k* dependent λ parameters
    TString klambda_filename = filename;
    auto file = std::make_unique<TFile>(klambda_filename, "read");
    auto lambdas = std::make_unique<TDirectory*[]>(5);
    lambdas[0] = (TDirectory*) file->GetDirectory("genuine");
    lambdas[1] = (TDirectory*) file->GetDirectory("FD_lambda");
    lambdas[2] = (TDirectory*) file->GetDirectory("FD_sigma");
    lambdas[3] = (TDirectory*) file->GetDirectory("flat");
    lambdas[4] = (TDirectory*) file->GetDirectory("fake");

    double scale = 0.1;
    for (size_t nlam = 0; nlam < 5; ++nlam)
    {
	TString name = Form("Mult_%i/Mt_%i/Lambda", var->mult, var->mt - 1);
	auto th1_klambda = (TH1F*) lambdas[nlam]->Get(name);
	th1_klambda->SetDirectory(0);

	if	(var->lam == 1) th1_klambda->Scale(1. - scale);
	else if (var->lam == 2) th1_klambda->Scale(1. + scale);

	lambda_pars[nlam] = th1_klambda;
    }
    file->Close();
}

void get_input_histos(TString ifile, TH1F **input_cf, TH1F **input_me, TH1F **input_me_orig)
{
    auto input_file = std::make_unique<TFile>(ifile, "read");
    *input_cf = (TH1F*) input_file->Get("CF");
    *input_me = (TH1F*) input_file->Get("ME");
    *input_me_orig = (TH1F*) input_file->Get("ME_original");
    (*input_cf)->SetDirectory(0);
    (*input_cf)->GetXaxis()->SetRangeUser(0., 700.);
    (*input_me)->SetDirectory(0);
    (*input_me_orig)->SetDirectory(0);
    input_file->Close();
}

void get_input_histos(TString ifile, std::unique_ptr<TH1F> &input_cf,
	std::unique_ptr<TH1F> &input_me, std::unique_ptr<TH1F> &input_me_orig)
{
    auto input_file = std::make_unique<TFile>(ifile, "read");
    input_cf.reset(static_cast<TH1F*>(input_file->Get("CF")));
    input_me.reset(static_cast<TH1F*>(input_file->Get("ME")));
    input_me_orig.reset(static_cast<TH1F*>(input_file->Get("ME_original")));
    input_cf->SetDirectory(0);
    input_cf->GetXaxis()->SetRangeUser(0., 700.);
    input_me->SetDirectory(0);
    input_me_orig->SetDirectory(0);
    input_file->Close();
}

void get_input_pl(TString ifile, TH1F **input_cf, TH1F **input_me, TH1F **input_me_orig, VAR *var, int charge)
{
    auto input_file = std::make_unique<TFile>(ifile, "read");

    TString mt[] = {"1.02-1.14", "1.14-1.2", "1.2-1.26", "1.26-1.38", "1.38-1.56", "1.56-1.86", "1.86-6"};
    TString mult[] = {"0.0-10.0", "10.0-50.0", "50.0-100.0"};

    TString tdir_name = (charge == APAP) ? "APAL/" : "PL/";
    tdir_name += "Rebin_8_Dim_1-" + mt[var->mt] + "_Dim_3-" + mult[var->mult];

    std::unique_ptr<TDirectory> tdir (static_cast<TDirectory*>(input_file->Get(tdir_name)));
    *input_cf = static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("CF/CF_rescaled"))->Clone("CF"));
    *input_me = static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("ME/ME_Reweighted"))->Clone("ME"));
    *input_me_orig = static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("No_Rebin/ME/ME"))->Clone("ME_original"));

    (*input_cf)->GetXaxis()->SetRangeUser(0., 700.);

    (*input_cf)->SetDirectory(0);
    (*input_me)->SetDirectory(0);
    (*input_me_orig)->SetDirectory(0);
}

void get_input_pl(TString ifile, std::unique_ptr<TH1F> &input_cf,
	std::unique_ptr<TH1F> &input_me, std::unique_ptr<TH1F> &input_me_orig, VAR *var, int charge)
{
    auto input_file = std::make_unique<TFile>(ifile, "read");

    TString mt[] = {"1.02-1.14", "1.14-1.2", "1.2-1.26", "1.26-1.38", "1.38-1.56", "1.56-1.86", "1.86-6"};
    TString mult[] = {"0.0-10.0", "10.0-50.0", "50.0-100.0"};

    TString tdir_name = (charge == APAP) ? "APAL/" : "PL/";
    tdir_name += "Rebin_8_Dim_1-" + mt[var->mt] + "_Dim_3-" + mult[var->mult];

    std::unique_ptr<TDirectory> tdir (static_cast<TDirectory*>(input_file->Get(tdir_name)));
    input_cf.reset(static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("CF/CF_rescaled"))->Clone("CF")));
    input_me.reset(static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("ME/ME_Reweighted"))->Clone("ME")));
    input_me_orig.reset(static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("No_Rebin/ME/ME"))->Clone("ME_original")));

    input_cf->GetXaxis()->SetRangeUser(0., 700.);

    input_cf->SetDirectory(0);
    input_me->SetDirectory(0);
    input_me_orig->SetDirectory(0);
}

void get_sample_histos(TString fname_cf, TString fname_syst, TH1F **input_cf, TH1F **input_me, TH1F **input_me_orig, VAR *var, int charge)
{
    get_input_histos(fname_cf, input_cf, input_me, input_me_orig);
    auto file_syst = std::make_unique<TFile>(fname_syst, "read");
    TString name_syst = Form("femto-dream-pair-task-track-track_std/mt_%i/mult_1/rebin_8/syst_th1", var->mt);

    std::unique_ptr<TH1F> input_syst (static_cast<TH1F*>(file_syst->Get(name_syst)));
    std::unique_ptr<TH1F> sample_cf (static_cast<TH1F*>((*input_cf)->Clone(Form("CF_sample_%s", CHARGE_STR[charge].Data()))));

    sample_cf->SetTitle(sample_cf->GetName());
    sample_cf->Reset();

    TRandom *random = new TRandom3(var->sample);
    double gaus, uniform, mean, sigma, syst;
    for (size_t n = 1; n < sample_cf->GetNbinsX(); ++n)
    {
	mean  = (*input_cf)->GetBinContent(n);
	sigma = (*input_cf)->GetBinError(n);
	syst  = input_syst->GetBinContent(n);

	gaus = random->Gaus(mean, sigma);
	uniform = random->Uniform(gaus - (syst / 2), gaus + (syst / 2));

	sample_cf->SetBinContent(n, uniform);
	if (var->stat)
	    sample_cf->SetBinContent(n, gaus);
	sample_cf->SetBinError(n, sigma);
    }
    *input_cf = static_cast<TH1F*>(sample_cf->Clone());
    (*input_cf)->SetDirectory(0);
}

void get_sample_histos(TString fname_cf, TString fname_syst, std::unique_ptr<TH1F> &input_cf,
	std::unique_ptr<TH1F> &input_me, std::unique_ptr<TH1F> &input_me_orig, VAR *var, int charge)
{
    get_input_histos(fname_cf, input_cf, input_me, input_me_orig);
    auto file_syst = std::make_unique<TFile>(fname_syst, "read");
    TString name_syst = Form("femto-dream-pair-task-track-track_std/mt_%i/mult_1/rebin_8/syst_th1", var->mt);

    std::unique_ptr<TH1F> input_syst (static_cast<TH1F*>(file_syst->Get(name_syst)));
    std::unique_ptr<TH1F> sample_cf (static_cast<TH1F*>(input_cf->Clone(Form("CF_sample_%s", CHARGE_STR[charge].Data()))));

    sample_cf->SetTitle(sample_cf->GetName());
    sample_cf->Reset();

    TRandom *random = new TRandom3(var->sample);
    double gaus, uniform, mean, sigma, syst;
    for (size_t n = 1; n < sample_cf->GetNbinsX(); ++n)
    {
	mean  = input_cf->GetBinContent(n);
	sigma = input_cf->GetBinError(n);
	syst  = input_syst->GetBinContent(n);

	gaus = random->Gaus(mean, sigma);
	uniform = random->Uniform(gaus - (syst / 2), gaus + (syst / 2));

	sample_cf->SetBinContent(n, uniform);
	if (var->stat)
	    sample_cf->SetBinContent(n, gaus);
	sample_cf->SetBinError(n, sigma);
    }
    input_cf.reset(static_cast<TH1F*>(sample_cf->Clone()));
    input_cf->SetDirectory(0);
}

template <typename T>
void get_input(TString input_path, TString fname_cf, T &input_cf, T &input_me, T &input_me_orig, VAR *var, int charge)
{
    if (var->sample)
    {
	TString name_syst = input_path + Form("UFFA_syst_%s_MultPercentile%i.root", CHARGE_STR[charge].Data(), var->mult);
	get_sample_histos(fname_cf, name_syst, input_cf, input_me, input_me_orig, var, charge);
    }
    else
    {
	get_input_histos(fname_cf, input_cf, input_me, input_me_orig);
    }
}

void get_input_histos_uffa(int nmt, int system, TH1F *&input_cf, TH1F *&input_me, TH1F *&input_me_orig)
{
    TString ifile = "UFFA_FullTrain_";
    ifile += (system == PP)? "0-100.root" : "apBase_0-100.root";
    auto input_file = std::make_unique<TFile>(ifile, "read");

    TString tdir_name = "femto-dream-pair-task-track-track_";
    tdir_name += (system == PP)? "id13610" : "apBase_id13610";
    tdir_name += Form("/bin_mt_%i/bin_1/");

    input_cf = (TH1F*) input_file->Get(tdir_name + "rebin_8/CF");
    input_me = (TH1F*) input_file->Get(tdir_name + "rebin_8/ME");
    input_me_orig = (TH1F*) ((TH1F*) input_file->Get(tdir_name + "ME"))->Clone("ME_original");
    input_cf->SetDirectory(0);
    input_cf->GetXaxis()->SetRangeUser(0., 700.);
    input_me->SetDirectory(0);
    input_me_orig->SetDirectory(0);
    input_file->Close();
}

void get_resolution_matrix(TString ipath, TH2F *&matrix, VAR *var, int charge)
{
    TString name_file, name_matrix;
    TFile *file;

    if (var->smear == 1)
    {
        name_file = "Resolution.root";
	name_matrix = "resolution_matrix_" + CHARGE_STR[charge];
    }
    else if (var->smear == 2)
    {
        name_file = "ALICE_pp_13TeV_MEpp.root";
        name_matrix = "h_RESO_pp_MeV";
    }

    if (var->smear)
    {
	auto file = std::make_unique<TFile>(ipath + "input/shared/" + name_file, "read");
	matrix = (TH2F*) file->Get(name_matrix);
        matrix->SetDirectory(0);

	printf("\e[1;34m  ┌───────────────┐\e[0m\n");
	printf("\e[1;34m  │  Res Matrix   │   \e[0m%s\n", name_file.Data());
	printf("\e[1;34m  │  Hist Name    │   \e[0m%s\n", name_matrix.Data());
	printf("\e[1;34m  └───────────────┘\e[0m\n\n");

	file->Close();
    }
    else
	matrix = NULL;
}

void get_cats_potential(CATS *cats, TH1F **pot_hist, TString target)
{
    double pot_bins = 1000, pot_min = 0, pot_max = 8;
    double radius;
    int channel = 0, part_wave = 0;

    if	    (!strcmp("1S0", target)) { channel = 0; part_wave = 0; }
    else if (!strcmp("1P1", target)) { channel = 0; part_wave = 1; }
    else if (!strcmp("3P0", target)) { channel = 1; part_wave = 1; }
    else if (!strcmp("3P1", target)) { channel = 2; part_wave = 1; }
    else if (!strcmp("3P2", target)) { channel = 3; part_wave = 1; }
    else if (!strcmp("1D2", target)) { channel = 0; part_wave = 2; }
    else    out("Partial Wave \"" << target << "\" not supported!");

    auto potential = std::make_unique<TH1F>("hPot", "hPot", pot_bins, pot_min, pot_max);
    for (size_t nbin = 1; nbin < pot_bins; ++nbin)
    {
	radius = potential->GetBinCenter(nbin);
	potential->SetBinContent(nbin, cats->EvaluateThePotential(channel, part_wave, 50, radius));
    }

    *pot_hist = (TH1F*) potential->Clone();
    (*pot_hist)->SetDirectory(0);
}

template <typename T>
void get_cats_potentials(CATS *cats, TH1F *pot_histos[], std::vector<T> pot_names)
{
    for (size_t npot = 0; npot < pot_names.size(); ++npot)
    {
	get_cats_potential(cats, &pot_histos[npot], pot_names[npot]);
	pot_histos[npot]->SetName(TString("hPot_") + pot_names[npot]);
	pot_histos[npot]->SetTitle(TString("hPot_") + pot_names[npot]);
    }
}

void setup_epos(DLM_CleverMcLevyResoTM *MagicSource, VAR_RSM *var_rsm, int restype, const char *file)
{
    auto epos_file = std::make_unique<TFile>(file, "read");
    std::unique_ptr<TNtuple> tntuple (static_cast<TNtuple*>(epos_file->Get("InfoTuple_ClosePairs")));

    Float_t k_D, fP1, fP2, fM1, fM2, Tau1, Tau2;
    Float_t AngleRcP1, AngleRcP2, AngleP1P2;
    double RanVal1, RanVal2, RanVal3;
    DLM_Random RanGen(11);

    std::vector<std::pair<std::string, Float_t *>> values = {std::pair{"k_D", &k_D},
	std::pair{"P1", &fP1}, std::pair{"P2", &fP2},
	std::pair{"M1", &fM1}, std::pair{"M2", &fM2},
	std::pair{"Tau1", &Tau1}, std::pair{"Tau2", &Tau2},
	std::pair{"AngleRcP1", &AngleRcP1}, std::pair{"AngleRcP2", &AngleRcP2}, std::pair{"AngleP1P2", &AngleP1P2}
    };

    for (auto branch : values)
    {
	tntuple->SetBranchAddress(branch.first.data(), branch.second);
    }

    for (size_t uEntry = 0; tntuple->GetEntry(uEntry); ++uEntry)
    {
	if (k_D > CUTOFF_Q) continue;

	RanVal1 = RanGen.Exponential(var_rsm->mass / (fP1 * var_rsm->tau));
	RanVal2 = RanGen.Exponential(var_rsm->mass2 / (fP2 * var_rsm->tau2));

	if	(restype == PR) MagicSource->AddBGT_PR(RanVal2, cos(AngleRcP2));
	else if (restype == RP) MagicSource->AddBGT_RP(RanVal1, -cos(AngleRcP1));
	else if (restype == RR) MagicSource->AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
    }
}

void setup_rsm_pp(DLM_CleverMcLevyResoTM *MagicSource, VAR_RSM *var_rsm)
{
    /* RSM Source Setup */
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, var_rsm->frac);
    MagicSource->SetUpReso(1, var_rsm->frac);
    MagicSource->InitNumMcIter(1000000);

    setup_epos(MagicSource, var_rsm, PR, "input/filesCATS/Source/EposAngularDist/EposDisto_p_pReso.root");
    setup_epos(MagicSource, var_rsm, RR, "input/filesCATS/Source/EposAngularDist/EposDisto_pReso_pReso.root");
}

void setup_rsm_pl(DLM_CleverMcLevyResoTM *MagicSource, VAR_RSM *var_rsm)
{
    /* RSM Source Setup */
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, var_rsm->frac);
    MagicSource->SetUpReso(1, var_rsm->frac2);
    MagicSource->InitNumMcIter(1000000);

    setup_epos(MagicSource, var_rsm, PR, "input/filesCATS/Source/EposAngularDist/EposDisto_p_LamReso.root");
    setup_epos(MagicSource, var_rsm, RP, "input/filesCATS/Source/EposAngularDist/EposDisto_pReso_Lam.root");
    setup_epos(MagicSource, var_rsm, RR, "input/filesCATS/Source/EposAngularDist/EposDisto_pReso_LamReso.root");
}

void setup_cats(DLM_CommonAnaFunctions *setupper, CATS *cats, VAR_FMR *range, const char *target, DLM_CleverMcLevyResoTM *MagicSource, VAR_RSM *var_rsm)
{
    cats->SetMomBins(range->Bins, range->Min, range->Max);
    if	    (!strcmp("pp", target))	{ setupper->SetUpCats_pp(*cats, POT_PP,   SOURCE, 0, 0); }
    else if (!strcmp("reid93", target)) { setupper->SetUpCats_pp(*cats, "ReidV8", SOURCE, 0, 0); }
    else if (!strcmp("reidSC", target)) { setupper->SetUpCats_pp(*cats, "ReidSC", SOURCE, 0, 0); }
    else if (!strcmp("psp", target))	{ setupper->SetUpCats_pSp(*cats, POT_PS, SOURCE, 0, 0); }
    else if (!strcmp("ps0", target))	{ setupper->SetUpCats_pS0(*cats, "Chiral", SOURCE); }
    else if (!strcmp("pxm", target))	{ setupper->SetUpCats_pXim(*cats, "pXim_HALQCDPaper2020", SOURCE); }
    else if (!strcmp("px0", target))	{ setupper->SetUpCats_pXi0(*cats, "pXim_HALQCDPaper2020", SOURCE); }
    else if (!strcmp("px1530", target)) { setupper->SetUpCats_pXim(*cats, "pXim1530", SOURCE); }
    else if (strstr(target, "av18"))
    {
	auto cPars = std::make_unique<CATSparameters>(CATSparameters::tSource, 1, true);
	cPars->SetParameter(0, 1.2);

	cats->SetMomentumDependentSource(false);
	cats->SetThetaDependentSource(false);
	cats->SetAnaSource(GaussSource, *cPars);
	cats->SetUseAnalyticSource(true);
	cats->SetExcludeFailedBins(false);

	// #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
	int POT_FLAG = v18_Coupled3P2;
	double PotPars1S0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 0, 0, 0};
	double PotPars3P0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 1, 1, 0};
	double PotPars3P1[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 1, 1, 1};
	double PotPars3P2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 1, 1, 2};
	double PotPars1D2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 0, 2, 2};

	auto cPotPars1S0 = std::make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);
	auto cPotPars3P0 = std::make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);
	auto cPotPars3P1 = std::make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);
	auto cPotPars3P2 = std::make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);
	auto cPotPars1D2 = std::make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);

	cPotPars1S0->SetParameters(PotPars1S0);
	cPotPars3P0->SetParameters(PotPars3P0);
	cPotPars3P1->SetParameters(PotPars3P1);
	cPotPars3P2->SetParameters(PotPars3P2);
	cPotPars1D2->SetParameters(PotPars1D2);

	cats->SetQ1Q2(1);
	cats->SetPdgId(2212, 2212);
	cats->SetRedMass(0.5 * Mass_p);

	if (!cats->GetNumChannels())
	{
	    cats->SetNumChannels(4);
	    cats->SetNumPW(0, 1);
	    cats->SetNumPW(1, 2);
	    cats->SetNumPW(2, 2);
	    cats->SetNumPW(3, 2);
	    cats->SetSpin(0, 0);
	    cats->SetSpin(1, 1);
	    cats->SetSpin(2, 1);
	    cats->SetSpin(3, 1);
	    cats->SetChannelWeight(0, 3. / 12.);
	    cats->SetChannelWeight(1, 1. / 12.);
	    cats->SetChannelWeight(2, 3. / 12.);
	    cats->SetChannelWeight(3, 5. / 12.);
	}

	if (strstr(target, "av18_s"))
	    cats->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
	if (strstr(target, "av18_sp"))
	{
	    cats->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
	    cats->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
	    cats->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
	}
    }
    else if (!strcmp("reid68", target))
    {
	CATSparameters *cats_pars = NULL;

	CATSparameters *cPotPars1S0 = NULL;
	CATSparameters *cPotPars1P1 = NULL;

	CATSparameters *cPotPars3S1 = NULL;
	CATSparameters *cPotPars3P0 = NULL;
	CATSparameters *cPotPars3P1 = NULL;
	CATSparameters *cPotPars3P2 = NULL;
	CATSparameters *cPotPars1D2 = NULL;

	cats_pars = new CATSparameters(CATSparameters::tSource, 1, true);
	cats_pars->SetParameter(0, 1.2);
	cats->SetAnaSource(GaussSource, *cats_pars);
	cats->SetUseAnalyticSource(true);

	double PotPars1S0[8] = {pp_ReidOli, 0};
	double PotPars3P0[8] = {pp_ReidOli, 1};
	double PotPars3P1[8] = {pp_ReidOli, 1};
	double PotPars3P2[8] = {pp_ReidOli, 1};
	cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 2, true);
	cPotPars1S0->SetParameters(PotPars1S0);
	cPotPars3P0 = new CATSparameters(CATSparameters::tPotential, 2, true);
	cPotPars3P0->SetParameters(PotPars3P0);
	//cPotPars3P1 = new CATSparameters(CATSparameters::tPotential, 2, true);
	//cPotPars3P1->SetParameters(PotPars3P1);
	//cPotPars3P2 = new CATSparameters(CATSparameters::tPotential, 2, true);
	//cPotPars3P2->SetParameters(PotPars3P2);

	cats->SetMomentumDependentSource(false);
	cats->SetExcludeFailedBins(false);

	cats->SetQ1Q2(1);
	cats->SetPdgId(2212, 2212);
	cats->SetRedMass(0.5 * Mass_p);

	cats->SetNumChannels(4);
	cats->SetNumPW(0, 3);
	cats->SetNumPW(1, 2);
	cats->SetNumPW(2, 2);
	cats->SetNumPW(3, 2);
	cats->SetSpin(0, 0);
	cats->SetSpin(1, 1);
	cats->SetSpin(2, 1);
	cats->SetSpin(3, 1);
	cats->SetChannelWeight(0, 3. / 12.);
	cats->SetChannelWeight(1, 1. / 12.);
	cats->SetChannelWeight(2, 3. / 12.);
	cats->SetChannelWeight(3, 5. / 12.);

	if (cPotPars1S0) cats->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
	if (cPotPars1P1) cats->SetShortRangePotential(3, 1, fDlmPot, *cPotPars1P1);
	if (cPotPars1D2) cats->SetShortRangePotential(0, 2, fDlmPot, *cPotPars1D2);
	if (cPotPars3S1) cats->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3S1);
	if (cPotPars3S1 && cPotPars3P1) cats->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3S1);
	if (cPotPars3S1 && cPotPars3P2) cats->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
	if (cPotPars3P0) cats->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
	if (cPotPars3P1) cats->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
	if (cPotPars3P2) cats->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
    }
    else if (!strcmp("rsm", target))
    {
	setup_rsm_pp(MagicSource, var_rsm);
	cats->SetAnaSource(CatsSourceForwarder, MagicSource, 2);
	cats->SetAnaSource(0, 1.0);
	cats->SetAnaSource(1, 2.0);
	cats->SetUseAnalyticSource(true);
	setupper->SetUpCats_pp(*cats, POT_PP, "", 0, 0);
    }
    else if (!strcmp("pl", target))
    {
	setupper->SetUpCats_pL(*cats, POT_PL, SOURCE, 11600, 0);
	cats->SetChannelWeight( 7,  1./4. * CUSP_WEIGHT);	//1S0 SN(s) -> LN(s)
	cats->SetChannelWeight( 8,  3./4. * CUSP_WEIGHT);	//3S1 SN(s) -> LN(s)
	cats->SetChannelWeight(10,  3./4. * CUSP_WEIGHT);	//3S1 SN(d) -> LN(s)
	cats->SetChannelWeight(13, 3./20. * CUSP_WEIGHT);	//3D1 SN(d) -> LN(d)
	cats->SetChannelWeight(15, 3./20. * CUSP_WEIGHT);	//3D1 SN(s) -> LN(d)
    }
    else if (!strcmp("plrsm", target))
    {
	setup_rsm_pl(MagicSource, var_rsm);
	cats->SetAnaSource(CatsSourceForwarder, MagicSource, 2);
	cats->SetAnaSource(0, 1.0);
	cats->SetAnaSource(1, 2.0);
	cats->SetUseAnalyticSource(true);

	setupper->SetUpCats_pL(*cats, POT_PL, "", 11600, 0);
	cats->SetChannelWeight( 7,  1./4. * CUSP_WEIGHT);	//1S0 SN(s) -> LN(s)
	cats->SetChannelWeight( 8,  3./4. * CUSP_WEIGHT);	//3S1 SN(s) -> LN(s)
	cats->SetChannelWeight(10,  3./4. * CUSP_WEIGHT);	//3S1 SN(d) -> LN(s)
	cats->SetChannelWeight(13, 3./20. * CUSP_WEIGHT);	//3D1 SN(d) -> LN(d)
	cats->SetChannelWeight(15, 3./20. * CUSP_WEIGHT);	//3D1 SN(s) -> LN(d)
    }

    cats->KillTheCat();
}

void setup_cats(DLM_CommonAnaFunctions *setupper, CATS *cats, VAR_FMR *range, const char *target)
{
    setup_cats(setupper, cats, range, target, NULL, NULL);
}

void setup_decomp(DLM_CkDecomposition *&decomp, DLM_Ck *&ck, CATS *cats, TH2F *res_matrix, VAR_FMR *range, const char *target)
{
    ck = new DLM_Ck(1, 0, *cats, range->Bins*2, range->Min, range->Max*2);
    ck->SetCutOff(range->Max, 700);
    ck->Update();

    const char *var1;
    int var2;
    if	    (!strcmp("pp",  target)) { var1 = "pp_Fitter";  var2 = 4; }
    else if (!strcmp("ppl", target)) { var1 = "pLambda";    var2 = 1; }
    else if (!strcmp("pl",  target)) { var1 = "pLambda";    var2 = 5; }
    else if (!strcmp("psp", target)) { var1 = "pSigmaPlus"; var2 = 1; }
    else if (!strcmp("ps0", target)) { var1 = "pSigma0";    var2 = 1; }
    else if (!strcmp("px0", target)) { var1 = "pXi0";	    var2 = 2; }
    else if (!strcmp("pxm", target)) { var1 = "pXim";	    var2 = 2; }
    else if (!strcmp("px1530", target)) { var1 = "pXim1530", var2 = 0; }

    decomp = new DLM_CkDecomposition(var1, var2, *ck, res_matrix);
}

void setup_fitter(TF1 *fitter, VAR *var, double radius)
{
    fitter->SetParameter(0, 0.98);
    fitter->FixParameter(1, 0.);
    fitter->SetParameter(2, 0.);
    fitter->SetParameter(3, 0.);
    fitter->FixParameter(4, 0.);

    if (var->bsl) fitter->FixParameter(3, 0.);

    fitter->SetParameter(5, radius);
    fitter->SetParLimits(5, 0.4, 2);
}

void setup_global_fitter(ROOT::Fit::Fitter *fitter, VAR *var, double global_pars[])
{
    fitter->Config().SetParamsSettings(11, global_pars);
    fitter->Config().ParSettings(0).SetName("N pp");
    fitter->Config().ParSettings(1).SetName("a pp");
    fitter->Config().ParSettings(2).SetName("b pp");
    fitter->Config().ParSettings(3).SetName("c pp");
    fitter->Config().ParSettings(4).SetName("d pp");
    fitter->Config().ParSettings(5).SetName("N apap");
    fitter->Config().ParSettings(6).SetName("a apap");
    fitter->Config().ParSettings(7).SetName("b apap");
    fitter->Config().ParSettings(8).SetName("c apap");
    fitter->Config().ParSettings(9).SetName("d apap");
    fitter->Config().ParSettings(10).SetName("Radius");
    double radius = INIT_RADIUS[var->rsm][var->mult][var->mt];
    fitter->Config().ParSettings(10).SetValue(radius);
    //fitter->Config().ParSettings(10).Fix();

    fitter->Config().ParSettings(1).Fix();
    fitter->Config().ParSettings(4).Fix();
    fitter->Config().ParSettings(6).Fix();
    fitter->Config().ParSettings(9).Fix();

    fitter->Config().ParSettings(10).SetLimits(0.4, 2.0);

    for (size_t n = 1; n < 10; ++n)
    {
	if (n == 5) continue;
	fitter->Config().ParSettings(n).SetLimits(-1., 1.);
    }

    if (var->bsl)
    {
	fitter->Config().ParSettings(3).SetValue(0);
	fitter->Config().ParSettings(8).SetValue(0);
	fitter->Config().ParSettings(3).Fix();
	fitter->Config().ParSettings(8).Fix();
    }

    fitter->Config().MinimizerOptions().SetPrintLevel(0);
    fitter->Config().SetMinimizer("Minuit2", "Migrad");
}

void create_output_2(TFile *&file, TString iname, TString opath, VAR *var)
{
    // Prepend "FitResults_" to inputfile name
    TObjArray* split_fname = iname.Tokenize("/");
    int split_length = split_fname->GetEntries();
    TString fname = "FitResults_" + ((TObjString*) (split_fname->At(split_length - 1)))->String();

    // Append variations to inputfile name
    split_fname = fname.Tokenize(".");
    TString variation = Form("_smearing_%i_norm_%i_femtorange_%i_fitrange_%i_lambda_%i_bl_%i_rsm_%i.root",
	    var->smear, var->normal, var->fmr, var->fr, var->lam, var->bsl, var->rsm);
    TString outputfile = ((TObjString*)(split_fname->At(0)))->String() + variation;

    file = new TFile(opath + outputfile, "recreate");
    printf("\n  \e[1;36m-->  \e[0;36mFile Created  \e[1;36m<--\e[0m\n\n");
    std::cout << "  " << opath + outputfile << "\n\n";
}

void create_output_default(TFile *&file, TString iname, TString opath, VAR *var)
{
    TString fname = iname + CHARGE_STR[var->charge];
    fname += (var->rsm)? "_rsm" : NULL;
    fname += Form("_mt%i_mult%i_fmr%i_fr%i_lam%i_bsl%i_smear%i",
	    var->mt, var->mult, var->fmr, var->fr, var->lam, var->bsl, var->smear);
    fname += (var->rsm)? Form("_frac%i_mass%i.root", var->frac, var->mass) : ".root";

    file = new TFile(opath + fname, "recreate");
    printf("  \e[1;36m-->  \e[0;36mFile Created  \e[1;36m<--\e[0m\n\n");
    std::cout << "  " << opath << "\e[1;34m" << fname << "\e[0m\n\n";
}

void create_output_sample(TFile *&file, TString iname, TString opath, VAR *var)
{
    TString fname = iname + CHARGE_STR[var->charge];
    fname += (var->rsm)? "_rsm" : NULL;
    if (var->sample && var->stat)
	fname += "_stat";
    fname += Form("_mt%i_mult%i_seed%i.root", var->mt, var->mult, var->sample);

    if (var->charge == 0) opath += "pp/";
    if (var->charge == 1) opath += "apap/";
    if (var->stat) opath += "stat/";
    if (var->rsm)  opath += "rsm/";

    file = new TFile(opath + fname, "recreate");
    printf("  \e[1;36m-->  \e[0;36mFile Created  \e[1;36m<--\e[0m\n\n");
    std::cout << "  " << opath << "\e[1;34m" << fname << "\e[0m\n\n";
}

void create_output(TFile *&file, TString iname, TString opath, VAR *var)
{
    if (var->sample)
	create_output_sample(file, iname, opath, var);
    else
	create_output_default(file, iname, opath, var);
}

void get_corrected_cf(TH1F *cfs[], TH1F *lambdas[], double bsl_pars[], TH1F **corrected_cf)
{
    enum CF_ID { RAW = 0, LAM = 1, SIG = 2 };
    enum LAM_ID { GEN = 0, FLAT = 3, FAKE = 4 };

    double value_x, value_y, bsl, raw, lam, sig, flat;
    double lam_gen, lam_lam, lam_sig, lam_flat, lam_fake;
    for (size_t nbin = 1; nbin < (*corrected_cf)->GetNbinsX(); ++nbin)
    {
	value_x = cfs[RAW]->GetBinCenter(nbin);
	bsl = bsl_pars[0] * (1 + bsl_pars[2]*pow(value_x, 2) + bsl_pars[3]*pow(value_x, 3));

	raw = cfs[RAW]->GetBinContent(nbin);
	lam = cfs[LAM]->GetBinContent(nbin);
	sig = cfs[SIG]->GetBinContent(nbin);
	lam_gen = lambdas[GEN]->GetBinContent(nbin);
	lam_lam = lambdas[LAM]->GetBinContent(nbin);
	lam_sig = lambdas[SIG]->GetBinContent(nbin);
	lam_flat = lambdas[FLAT]->GetBinContent(nbin);
	lam_fake = lambdas[FAKE]->GetBinContent(nbin);

	value_y = (1/lam_gen)*((raw/bsl) - lam_lam*lam - lam_sig*sig - lam_flat - lam_fake);

	if (value_y > 10e-6 && value_y < 50)
	    (*corrected_cf)->SetBinContent(nbin, value_y);
    }
}

void cf_fitter(TString project_path, TString file_name, TString output_path, VAR *var)
{
    printf("\n  \e[1;36m-->  \e[0;36mEntering the fitter!  \e[1;36m<--\e[0m\n\n");

    VAR_RSM *var_rsm = new VAR_RSM(PROTON, PROTON, var->rsm);
    VAR_FMR *range_femto = new VAR_FMR(var->system, var->fmr);
    VAR_FR  *range_fit = new VAR_FR(var->fr);
    VAR_QA  *qa_plots = new VAR_QA(var);

    if (var->mt > 6)
    {
	range_femto->Min = 8;
	range_fit->Min = 8;
    }
    if (var->mt == 4 && var->mult == 2)
    {
	range_femto->Min = 8;
	range_fit->Min = 8;
    }

    print_info_2(var, range_femto, range_fit);

    std::unique_ptr<TH1F> input_cf, input_me, input_me_orig;
    get_input_histos(file_name, input_cf, input_me, input_me_orig);

    std::unique_ptr<TH1F*> uptr_lambda_pars_th1 (new TH1F*[5]);
    auto lambda_pars_th1 = static_cast<TH1F**>(uptr_lambda_pars_th1.get());

    std::unique_ptr<DLM_Histo<double>*> uptr_lambda_pars (new DLM_Histo<double>*[5]);
    auto lambda_pars = static_cast<DLM_Histo<double>**>(uptr_lambda_pars.get());

    TString file_klambda = project_path + "input/shared/Lambda_" + ((var->charge)? "Anti" : "") + "Proton.root";
    get_klambda_histos(file_klambda, var, lambda_pars_th1);
    get_klambda_histos(file_klambda, var, lambda_pars);

    TH2F *matrix_pp_res;
    get_resolution_matrix(project_path, matrix_pp_res, var, var->charge);

    DLM_CommonAnaFunctions cats_setupper;
    cats_setupper.SetCatsFilesFolder(project_path + "input/filesCATS");
    TH2F *matrix_pl_res = cats_setupper.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
    TH2F *matrix_pl_dec = cats_setupper.GetResidualMatrix("pp", "pLambda");
    TH2F *matrix_ps_dec = cats_setupper.GetResidualMatrix("pp", "pSigmaPlus");

    DLM_CleverMcLevyResoTM MagicSource;
    DLM_CleverMcLevyResoTM *pMS = (var->rsm)? &MagicSource : NULL;
    VAR_RSM *pRSM = (var->rsm)? var_rsm : NULL;
    TString target = (var->rsm)? "rsm" : "pp";

    CATS cats_pp, cats_pl, cats_ps;
    setup_cats(&cats_setupper, &cats_pp, range_femto, target.Data(), pMS, pRSM);
    setup_cats(&cats_setupper, &cats_pl, range_femto, "pl");
    setup_cats(&cats_setupper, &cats_ps, range_femto, "psp");

    DLM_Ck *ck_pp, *ck_pl, *ck_ps;
    DLM_CkDecomposition *decomp_pp, *decomp_pl, *decomp_ps;
    setup_decomp(decomp_pp, ck_pp, &cats_pp, matrix_pp_res, range_femto, "pp");
    setup_decomp(decomp_pl, ck_pl, &cats_pl, matrix_pl_res, range_femto, "ppl");
    setup_decomp(decomp_ps, ck_ps, &cats_ps, matrix_pl_res, range_femto, "psp");

    // non-genuine contributions to pp
    decomp_pp->AddContribution(0, *lambda_pars[1], DLM_CkDecomp::cFeedDown, decomp_pl, matrix_pl_dec);
    decomp_pp->AddContribution(1, *lambda_pars[2], DLM_CkDecomp::cFeedDown, decomp_ps, matrix_ps_dec);
    decomp_pp->AddContribution(2, *lambda_pars[3], DLM_CkDecomp::cFeedDown);
    decomp_pp->AddContribution(3, *lambda_pars[4], DLM_CkDecomp::cFake);
    decomp_pp->AddPhaseSpace(input_me.get());
    decomp_pp->AddPhaseSpace(0, input_me.get());

    // non-genuine contributions to pL
    decomp_pl->AddContribution(0, LAM_PL_FLAT, DLM_CkDecomp::cFeedDown);

    // update CkDec objects
    decomp_pp->Update(true, true);
    decomp_pl->Update(true, true);
    decomp_ps->Update(true, true);

    FitFunSimpleAvg fitfun(input_cf.get(), input_me_orig.get(), decomp_pp);
    TF1 *fitter = new TF1("Fit", fitfun, range_fit->Min, range_fit->Max, 6);
    TF1 *baseline = new TF1("Baseline", bslfun, range_fit->Min, range_fit->Max, 5);

    setup_fitter(fitter, var, 1.);
    input_cf->Fit(fitter, "S, N, R, M");

    double bsl_pars[5];
    for (int npar = 0; npar < 5; ++npar)
    {
	bsl_pars[npar] = fitter->GetParameter(npar);
	baseline->SetParameter(npar, fitter->GetParameter(npar));
    }

    GenuineFun genfun(decomp_pp);
    TF1* genuine = new TF1("Genuine", genfun, range_fit->Min, range_fit->Max, 1);
    genuine->SetParameter(0, 1.);

    GenuineFunSmear genfun_smear(decomp_pp);
    TF1* genuine_smear = new TF1("Genuine Smeared", genfun_smear, range_fit->Min, range_fit->Max, 1);
    genuine_smear->SetParameter(0, 1.);

    decomp_pp->Update(true, true);
    TH1F *child_pl_smear = new TH1F("pl_smeared",   "pL smeared pp",    500, 0, 500);
    TH1F *child_ps_smear = new TH1F("ps_smeared",   "pS+ smeared pp",   500, 0, 500);
    for (short nbin = 0; nbin < 500; ++nbin)
    {
	child_pl_smear->SetBinContent(nbin + 1, decomp_pp->EvalSignalSmearedChild(0, child_pl_smear->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear->SetBinContent(nbin + 1, decomp_pp->EvalSignalSmearedChild(1, child_ps_smear->GetBinCenter(nbin + 1) + 1) + 1);
    }

    TH1F **cfs_corr = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_corr[0] = static_cast<TH1F*>(input_cf->Clone());
    cfs_corr[1] = child_pl_smear;
    cfs_corr[2] = child_ps_smear;

    TH1F *corrected_cf = static_cast<TH1F*>(input_cf->Clone("CF_corrected"));
    corrected_cf->Reset();
    get_corrected_cf(cfs_corr, lambda_pars_th1, bsl_pars, &corrected_cf);

    printf("\n\e[1;36m  -->  \e[0;36mData Fit  \e[1;36m<--\e[0m\n\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Source Size  │   \e[0m%.2f ± %.2f\n", fitter->GetParameter(5), fitter->GetParError(5));
    printf("\e[1;34m  │  p-λ Size     │   \e[0m%.2f\n", decomp_pp->GetChild(0)->GetCk()->GetSourcePar(0));
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");

    TFile *ofile;
    create_output(ofile, "FitResults_", output_path, var);
    ofile->cd();

    TString str_charge = CHARGE_STR[var->charge];
    input_cf->Write("CF_" + str_charge);
    corrected_cf->Write("CF_corr_" + str_charge);
    input_me_orig->Write("ME_orig_" + str_charge);

    fitter->Write("fFitResult_" + str_charge);
    baseline->Write("fBl_" + str_charge);

    genuine->Write("fGenuine_" + str_charge);
    genuine_smear->Write("fGenuineSmeared_" + str_charge);

    child_pl_smear->Write("pl_smeared_" + str_charge);
    child_ps_smear->Write("ps_smeared_" + str_charge);

    TDirectory *tdir_lam = ofile->mkdir("lam_pars");
    tdir_lam->cd();

    lambda_pars_th1[0]->Write("genuine_" + str_charge);
    lambda_pars_th1[1]->Write("lambda_" + str_charge);
    lambda_pars_th1[2]->Write("sigma_" + str_charge);
    lambda_pars_th1[3]->Write("flat_" + str_charge);
    lambda_pars_th1[4]->Write("fake_" + str_charge);

    TDirectory *tdir_qa = ofile->mkdir("qa_plots");
    tdir_qa->cd();

    qa_plots->lam->Write();
    qa_plots->fmr->Write();
    qa_plots->fr->Write();
    qa_plots->bsl->Write();
    qa_plots->smear->Write();

    if (var->rsm)
    {
	qa_plots->frac->Write();
	qa_plots->mass->Write();
    }

    ofile->Close();
}

void cf_fitter_pl(TString project_path, TString file_name, TString output_path, VAR *var)
{
    printf("\n  \e[1;36m-->  \e[0;36mEntering the fitter!  \e[1;36m<--\e[0m\n\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Input File   │\e[0m   %s\e[0m\n", file_name.Data());
    printf("\e[1;34m  └───────────────┘\e[0m\n");

    VAR_RSM *var_rsm = new VAR_RSM(PROTON, LAMBDA, var->rsm);
    VAR_FMR *range_femto = new VAR_FMR(var->system, var->fmr);
    VAR_FR  *range_fit = new VAR_FR(var->fr);
    VAR_QA  *qa_plots = new VAR_QA(var);

    range_femto->Max = 340;
    range_femto->Min = 0;

    range_fit->Max = 500;
    range_fit->Min = 0;

    print_info_2(var, range_femto, range_fit);

    TH1F *input_cf, *input_me, *input_me_orig;
    get_input_pl(file_name, &input_cf, &input_me, &input_me_orig, var, var->charge);

    /* CATS Setup */
    DLM_CommonAnaFunctions cats_setupper;
    cats_setupper.SetCatsFilesFolder(project_path + "input/filesCATS");
    TH2F *matrix_pl_res = cats_setupper.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
    TH2F *matrix_pl_dec = cats_setupper.GetResidualMatrix("pp", "pLambda");
    TH2F *matrix_ls0_dec = cats_setupper.GetResidualMatrix("pLambda", "pSigma0");
    TH2F *matrix_lxm_dec = cats_setupper.GetResidualMatrix("pLambda", "pXim");
    TH2F *matrix_lx1530_dec = cats_setupper.GetResidualMatrix("pXim", "pXim1530");

    //std::vector<float> cusp_var {0.40, 0.27, 0.33}; // variations of the cusp

    /* CATS Objects */
    CATS cats_pl, cats_ps0, cats_px0, cats_pxm, cats_px1530;
    setup_cats(&cats_setupper, &cats_pl,  range_femto, "pl");
    setup_cats(&cats_setupper, &cats_ps0, range_femto, "ps0");
    setup_cats(&cats_setupper, &cats_px0, range_femto, "px0");
    setup_cats(&cats_setupper, &cats_pxm, range_femto, "pxm");
    setup_cats(&cats_setupper, &cats_px1530, range_femto, "px1530");

    /* Decomp Objects */
    DLM_Ck *ck_pl, *ck_ps0, *ck_px0, *ck_pxm, *ck_px1530;
    DLM_CkDecomposition *decomp_pl, *decomp_ps0, *decomp_px0, *decomp_pxm, *decomp_px1530;
    setup_decomp(decomp_pl,  ck_pl,  &cats_pl,  matrix_pl_res,	range_femto, "pl");
    setup_decomp(decomp_ps0, ck_ps0, &cats_ps0,	NULL,		range_femto, "ps0");
    setup_decomp(decomp_px0, ck_px0, &cats_px0, NULL,		range_femto, "px0");
    setup_decomp(decomp_pxm, ck_pxm, &cats_pxm, NULL,		range_femto, "pxm");
    setup_decomp(decomp_px1530, ck_px1530, &cats_px1530, NULL,	range_femto, "px1530");

    std::unique_ptr<double[]> lam_pars_pl  (new double[5]);
    std::unique_ptr<double[]> lam_pars_pxi (new double[5]);
    cats_setupper.SetUpLambdaPars_pL("pp13TeV_HM_Dec19", 0, 0, lam_pars_pl.get());
    cats_setupper.SetUpLambdaPars_pXim("pp13TeV_HM_Dec19", 0, 0, lam_pars_pxi.get());

    double lam_pXi_gen = lam_pars_pxi[0] / (1. - lam_pars_pxi[4]);
    double lam_pXi_pXi1530 = 0.5*(lam_pars_pxi[1]+lam_pars_pxi[2])/(1.-lam_pars_pxi[4]);
    double lam_pXi_flt = 1. - lam_pXi_gen - lam_pXi_pXi1530;

    double lam_gen = lam_pars_pl[0];
    double lam_pS0 = lam_pars_pl[1];
    double lam_pXim = lam_pars_pl[2];
    double lam_pXi0 = 0;
    double lam_flt = lam_pars_pl[3];
    double lam_mid = lam_pars_pl[4];

    double lam_pS0_flt = 0.18;

    printf(" lam_gen = %.2f%%\n",lam_gen*100.);
    printf(" lam_pS0 = %.2f%%\n",lam_pS0*100.);
    printf(" lam_pXim= %.2f%%\n",lam_pXim*100.);
    printf(" lam_pXi0= %.2f%%\n",lam_pXi0*100.);
    printf(" lam_flt = %.2f%%\n",lam_flt*100.);
    printf(" lam_mid = %.2f%%\n",lam_mid*100.);
    printf(" lam_pXi_gen = %.2f%%\n",lam_pXi_gen*100.);
    printf(" lam_pXi_pXi1530 = %.2f%%\n",lam_pXi_pXi1530*100.);
    printf(" lam_pXi_flt = %.2f%%\n",lam_pXi_flt*100.);
    printf(" lam_pS0_flt = %.2f%%\n",lam_pS0_flt*100.);

    decomp_pl->AddContribution(0, lam_pS0,  DLM_CkDecomp::cFeedDown, decomp_ps0, matrix_ls0_dec);
    decomp_pl->AddContribution(1, lam_pXim, DLM_CkDecomp::cFeedDown, decomp_pxm, matrix_lxm_dec);
    decomp_pl->AddContribution(2, lam_pXi0, DLM_CkDecomp::cFeedDown, decomp_px1530, matrix_lx1530_dec);
    decomp_pl->AddContribution(3, lam_flt,  DLM_CkDecomp::cFeedDown);
    decomp_pl->AddContribution(4, lam_mid,  DLM_CkDecomp::cFake);

    decomp_ps0->AddContribution(1, lam_pS0_flt, DLM_CkDecomposition::cFeedDown);

    decomp_pxm->AddContribution(0, lam_pXi_pXi1530, DLM_CkDecomposition::cFeedDown, decomp_px1530, matrix_lx1530_dec);
    decomp_pxm->AddContribution(1, lam_pXi_flt, DLM_CkDecomposition::cFeedDown);

    decomp_px0->AddContribution(0, lam_pXi_pXi1530, DLM_CkDecomposition::cFeedDown, decomp_px1530, matrix_lx1530_dec);
    decomp_px0->AddContribution(1, lam_pXi_flt, DLM_CkDecomposition::cFeedDown);

    ck_pl->Update();
    ck_ps0->Update();
    ck_px0->Update();
    ck_pxm->Update();
    ck_px1530->Update();

    //decomp_pl->AddContribution(0, LAM_PL["ps0"],  DLM_CkDecomp::cFeedDown, decomp_ps0, matrix_ls0_dec);
    //decomp_pl->AddContribution(1, LAM_PL["xim"],  DLM_CkDecomp::cFeedDown, decomp_pxm, matrix_lxm_dec);
    //decomp_pl->AddContribution(2, LAM_PL["flat"], DLM_CkDecomp::cFeedDown);
    //decomp_pl->AddContribution(3, LAM_PL["fake"], DLM_CkDecomp::cFake);
    //decomp_pl->AddPhaseSpace(input_me);
    //decomp_pl->Update(true, true);

    /* Fit Objects */
    plFitFunctionSimple fitfun(decomp_pl);
    TF1 *fitter = new TF1("Fit", fitfun, range_fit->Min, range_fit->Max, 6);
    TF1 *baseline = new TF1("Baseline", bslfun, range_fit->Min, range_fit->Max, 5);

    setup_fitter(fitter, var, 1.);
    input_cf->Fit(fitter, "S, N, R, M");

    double bsl_pars[5];
    for (int npar = 0; npar < 5; ++npar)
    {
	bsl_pars[npar] = fitter->GetParameter(npar);
	baseline->SetParameter(npar, fitter->GetParameter(npar));
    }

    GenuineFun genfun(decomp_pl);
    TF1* genuine = new TF1("Genuine", genfun, range_fit->Min, range_fit->Max, 1);
    genuine->SetParameter(0, 1.);

    GenuineFunSmear genfun_smear(decomp_pl);
    TF1* genuine_smear = new TF1("Genuine Smeared", genfun_smear, range_fit->Min, range_fit->Max, 1);
    genuine_smear->SetParameter(0, 1.);

    decomp_pl->Update(true, true);

    printf("\n\e[1;36m  -->  \e[0;36mData Fit  \e[1;36m<--\e[0m\n\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Source Size  │   \e[0m%.2f ± %.2f\n", fitter->GetParameter(5), fitter->GetParError(5));
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");

    TFile *ofile;
    create_output(ofile, "FitResults_", output_path, var);
    ofile->cd();

    TString str_charge = CHARGE_STR[var->charge];
    input_cf->Write("CF_" + str_charge);
    input_me_orig->Write("ME_orig_" + str_charge);

    fitter->Write("fFitResult_" + str_charge);
    baseline->Write("fBl_" + str_charge);

    genuine->Write("fGenuine_" + str_charge);
    genuine_smear->Write("fGenuineSmeared_" + str_charge);

    TDirectory *tdir_qa = ofile->mkdir("qa_plots");
    tdir_qa->cd();

    qa_plots->lam->Write();
    qa_plots->fmr->Write();
    qa_plots->fr->Write();
    qa_plots->bsl->Write();
    qa_plots->smear->Write();

    if (var->rsm)
    {
	qa_plots->frac->Write();
	qa_plots->mass->Write();
    }

    ofile->Close();
}

void cf_combined_fitter(TString project_path, TString input_path, TString output_path, TString file_name_pp, TString file_name_aa, VAR *var)
{
    printf("\n  \e[1;36m-->  \e[0;36mEntering the combined fitter!  \e[1;36m<--\e[0m\n\n");

    VAR_RSM *var_rsm = new VAR_RSM(PROTON, PROTON, var->rsm);
    VAR_FMR *range_femto_pp = new VAR_FMR(var->system, var->fmr);
    VAR_FMR *range_femto_aa = new VAR_FMR(var->system, var->fmr);
    VAR_FR  *range_fit_pp = new VAR_FR(var->fr);
    VAR_FR  *range_fit_aa = new VAR_FR(var->fr);
    VAR_QA  *qa_plots = new VAR_QA(var);

    // This setup is because we are missing the first cf bin for these mt and mult bins
    if (var->mt > 6)
    {
	range_femto_pp->Min = 8;
	range_femto_aa->Min = 8;
	range_fit_pp->Min = 8;
	range_fit_aa->Min = 8;
    }
    if (var->mt == 4 && var->mult == 2)
    {
	range_femto_aa->Min = 8;
	range_fit_aa->Min = 8;
    }

    print_info_2(var, range_femto_pp, range_fit_pp);

    // fitter parameter entry list, i.e. the last parameter is shared between both fitters
    int fit_pars_pp[6] = { 0, 1, 2, 3, 4, 10 };	    // N, a, b, c, d, r
    int fit_pars_aa[6] = { 5, 6, 7, 8, 9, 10 };
    double global_pars[11] = {0.98, 0, 0, 0, 0, 0.98, 0, 0, 0, 0, 1};	    // init parameters

    std::unique_ptr<TH1F> input_cf_pp, input_me_pp, input_me_orig_pp;
    std::unique_ptr<TH1F> input_cf_aa, input_me_aa, input_me_orig_aa;
    get_input(input_path, file_name_pp, input_cf_pp, input_me_pp, input_me_orig_pp, var, PP);
    get_input(input_path, file_name_aa, input_cf_aa, input_me_aa, input_me_orig_aa, var, APAP);

    std::unique_ptr<TH1F*> uptr_lambda_pars_pp_th1 (new TH1F*[5]);
    std::unique_ptr<TH1F*> uptr_lambda_pars_aa_th1 (new TH1F*[5]);
    auto lambda_pars_pp_th1 = static_cast<TH1F**>(uptr_lambda_pars_pp_th1.get());
    auto lambda_pars_aa_th1 = static_cast<TH1F**>(uptr_lambda_pars_aa_th1.get());
    get_klambda_histos(project_path + "input/shared/Lambda_Proton.root", var, lambda_pars_pp_th1);
    get_klambda_histos(project_path + "input/shared/Lambda_AntiProton.root", var, lambda_pars_aa_th1);

    std::unique_ptr<DLM_Histo<double>*> uptr_lambda_pars_pp (new DLM_Histo<double>*[5]);
    std::unique_ptr<DLM_Histo<double>*> uptr_lambda_pars_aa (new DLM_Histo<double>*[5]);
    auto lambda_pars_pp = static_cast<DLM_Histo<double>**>(uptr_lambda_pars_pp.get());
    auto lambda_pars_aa = static_cast<DLM_Histo<double>**>(uptr_lambda_pars_aa.get());
    get_klambda_histos(project_path + "input/shared/Lambda_Proton.root", var, lambda_pars_pp);
    get_klambda_histos(project_path + "input/shared/Lambda_AntiProton.root", var, lambda_pars_aa);

    TH2F *matrix_pp_res, *matrix_aa_res;
    get_resolution_matrix(project_path, matrix_pp_res, var, PP);
    get_resolution_matrix(project_path, matrix_aa_res, var, APAP);

    DLM_CommonAnaFunctions cats_setupper;
    cats_setupper.SetCatsFilesFolder(project_path + "input/filesCATS");
    TH2F *matrix_pl_res = cats_setupper.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
    TH2F *matrix_pl_dec = cats_setupper.GetResidualMatrix("pp", "pLambda");
    TH2F *matrix_ps_dec = cats_setupper.GetResidualMatrix("pp", "pSigmaPlus");

    DLM_CleverMcLevyResoTM MagicSource_pp, MagicSource_aa;
    DLM_CleverMcLevyResoTM *pMS_pp = (var->rsm)? &MagicSource_pp : NULL;
    DLM_CleverMcLevyResoTM *pMS_aa = (var->rsm)? &MagicSource_aa : NULL;
    VAR_RSM *pRSM = (var->rsm)? var_rsm : NULL;
    TString target = (var->rsm)? "rsm" : "pp";

    CATS cats_pp, cats_pl_pp, cats_ps_pp;
    setup_cats(&cats_setupper, &cats_pp,    range_femto_pp, target.Data(), pMS_pp, pRSM);
    setup_cats(&cats_setupper, &cats_pl_pp, range_femto_pp, "pl");
    setup_cats(&cats_setupper, &cats_ps_pp, range_femto_pp, "psp");

    CATS cats_aa, cats_pl_aa, cats_ps_aa;
    setup_cats(&cats_setupper, &cats_aa,    range_femto_aa, target.Data(), pMS_aa, pRSM);
    setup_cats(&cats_setupper, &cats_pl_aa, range_femto_aa, "pl");
    setup_cats(&cats_setupper, &cats_ps_aa, range_femto_aa, "psp");

    DLM_Ck *ck_pp, *ck_pl_pp, *ck_ps_pp;
    DLM_CkDecomposition *decomp_pp, *decomp_pl_pp, *decomp_ps_pp;
    setup_decomp(decomp_pp,    ck_pp,	 &cats_pp,    matrix_pp_res, range_femto_pp, "pp");
    setup_decomp(decomp_pl_pp, ck_pl_pp, &cats_pl_aa, matrix_pl_res, range_femto_pp, "ppl");
    setup_decomp(decomp_ps_pp, ck_ps_pp, &cats_ps_aa, matrix_pl_res, range_femto_pp, "psp");

    DLM_Ck *ck_aa, *ck_pl_aa, *ck_ps_aa;
    DLM_CkDecomposition *decomp_aa, *decomp_pl_aa, *decomp_ps_aa;
    setup_decomp(decomp_aa,    ck_aa,	 &cats_aa,    matrix_aa_res, range_femto_aa, "pp");
    setup_decomp(decomp_pl_aa, ck_pl_aa, &cats_pl_aa, matrix_pl_res, range_femto_aa, "ppl");
    setup_decomp(decomp_ps_aa, ck_ps_aa, &cats_ps_aa, matrix_pl_res, range_femto_aa, "psp");

    // non-genuine contributions to pp
    decomp_pp->AddContribution(0, *lambda_pars_pp[1], DLM_CkDecomp::cFeedDown, decomp_pl_pp, matrix_pl_dec);
    decomp_pp->AddContribution(1, *lambda_pars_pp[2], DLM_CkDecomp::cFeedDown, decomp_ps_pp, matrix_ps_dec);
    decomp_pp->AddContribution(2, *lambda_pars_pp[3], DLM_CkDecomp::cFeedDown);
    decomp_pp->AddContribution(3, *lambda_pars_pp[4], DLM_CkDecomp::cFake);
    decomp_pp->AddPhaseSpace(input_me_pp.get());
    decomp_pp->AddPhaseSpace(0, input_me_pp.get());

    // non-genuine contributions to apap
    decomp_aa->AddContribution(0, *lambda_pars_aa[1], DLM_CkDecomp::cFeedDown, decomp_pl_aa, matrix_pl_dec);
    decomp_aa->AddContribution(1, *lambda_pars_aa[2], DLM_CkDecomp::cFeedDown, decomp_ps_aa, matrix_ps_dec);
    decomp_aa->AddContribution(2, *lambda_pars_aa[3], DLM_CkDecomp::cFeedDown);
    decomp_aa->AddContribution(3, *lambda_pars_aa[4], DLM_CkDecomp::cFake);
    decomp_aa->AddPhaseSpace(input_me_aa.get());
    decomp_aa->AddPhaseSpace(0, input_me_aa.get());

    // non-genuine contributions to pL
    decomp_pl_pp->AddContribution(0, LAM_PL_FLAT, DLM_CkDecomp::cFeedDown);
    decomp_pl_aa->AddContribution(0, LAM_PL_FLAT, DLM_CkDecomp::cFeedDown);

    // update CkDec objects
    decomp_pp->Update(true, true);
    decomp_pl_pp->Update(true, true);
    decomp_ps_pp->Update(true, true);

    decomp_aa->Update(true, true);
    decomp_pl_aa->Update(true, true);
    decomp_ps_aa->Update(true, true);

    // Fit function
    FitFunSimpleAvg fitfun_pp(input_cf_pp.get(), input_me_orig_pp.get(), decomp_pp);
    FitFunSimpleAvg fitfun_aa(input_cf_aa.get(), input_me_orig_aa.get(), decomp_aa);
    TF1 *fitter_pp = new TF1("Fit pp",	 fitfun_pp, range_fit_pp->Min, range_fit_pp->Max, 6);
    TF1 *fitter_aa = new TF1("Fit apap", fitfun_aa, range_fit_aa->Min, range_fit_aa->Max, 6);

    setup_fitter(fitter_pp, var, INIT_RADIUS[var->rsm][var->mult][var->mt]);
    setup_fitter(fitter_aa, var, INIT_RADIUS[var->rsm][var->mult][var->mt]);

    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange opt_range_pp, opt_range_aa;
    opt_range_pp.SetRange(range_fit_pp->Min, range_fit_pp->Max);
    opt_range_aa.SetRange(range_fit_aa->Min, range_fit_aa->Max);

    ROOT::Fit::BinData wrap_data_pp(opt, opt_range_pp);
    ROOT::Fit::BinData wrap_data_aa(opt, opt_range_aa);

    ROOT::Math::WrappedMultiTF1 wrap_fit_pp(*fitter_pp, 1);
    ROOT::Math::WrappedMultiTF1 wrap_fit_aa(*fitter_aa, 1);

    ROOT::Fit::FillData(wrap_data_pp, input_cf_pp.get());
    ROOT::Fit::FillData(wrap_data_aa, input_cf_aa.get());

    ROOT::Fit::Chi2Function chi2_pp(wrap_data_pp, wrap_fit_pp);
    ROOT::Fit::Chi2Function chi2_aa(wrap_data_aa, wrap_fit_aa);

    GlobalChi2 global_chi2(chi2_pp, chi2_aa, fit_pars_pp, fit_pars_aa);

    ROOT::Fit::Fitter fitter;
    setup_global_fitter(&fitter, var, global_pars);

    fitter.FitFCN(11, global_chi2, nullptr, wrap_data_pp.Size() + wrap_data_aa.Size(), true);
    ROOT::Fit::FitResult result = fitter.Result();

    fitter_pp->SetFitResult(result, fit_pars_pp);
    fitter_aa->SetFitResult(result, fit_pars_aa);

    double temp, total_chi2 = 0.;
    for (size_t nbin = input_cf_pp->FindBin(50); input_cf_pp->GetBinCenter(nbin) < 140; ++nbin)
    {
	temp = (std::pow(input_cf_pp->GetBinContent(nbin) - fitter_pp->Eval(input_cf_pp->GetBinCenter(nbin)), 2)) / std::pow(input_cf_pp->GetBinError(nbin), 2);
	total_chi2 += temp;
    }
    //std::cout << "pp fit Chi2 = " << total_chi2 << std::endl;
    out("pp fit Chi2 = " << total_chi2);

    fitter_pp->SetRange(opt_range_pp().first, opt_range_pp().second);
    fitter_aa->SetRange(opt_range_aa().first, opt_range_aa().second);

    double bsl_pars_pp[5], bsl_pars_aa[5];
    TF1 *baseline_pp = new TF1("Baseline pp",   bslfun, range_fit_pp->Min, range_fit_pp->Max, 5);
    TF1 *baseline_aa = new TF1("Baseline apap", bslfun, range_fit_aa->Min, range_fit_aa->Max, 5);
    for (size_t npar = 0; npar < 5; ++npar)
    {
	bsl_pars_pp[npar] = fitter_pp->GetParameter(npar);
	bsl_pars_aa[npar] = fitter_aa->GetParameter(npar);
	baseline_pp->SetParameter(npar, fitter_pp->GetParameter(npar));
	baseline_aa->SetParameter(npar, fitter_aa->GetParameter(npar));
    }

    GenuineFun genfun_pp(decomp_pp), genfun_aa(decomp_aa);
    TF1* genuine_pp = new TF1("Genuine pp", genfun_pp, range_fit_pp->Min, range_fit_pp->Max, 1);
    TF1* genuine_aa = new TF1("Genuine apap", genfun_aa, range_fit_aa->Min, range_fit_aa->Max, 1);
    genuine_aa->SetParameter(0, 1.);
    genuine_pp->SetParameter(0, 1.);

    GenuineFunSmear genfun_smear_pp(decomp_pp), genfun_smear_aa(decomp_aa);
    TF1* genuine_smear_pp = new TF1("Genuine Smeared pp", genfun_smear_pp, range_fit_pp->Min, range_fit_pp->Max, 1);
    TF1* genuine_smear_aa = new TF1("Genuine Smeared apap", genfun_smear_aa, range_fit_aa->Min, range_fit_aa->Max, 1);
    genuine_smear_pp->SetParameter(0, 1.);
    genuine_smear_aa->SetParameter(0, 1.);

    decomp_pp->Update(true, true);
    decomp_aa->Update(true, true);
    TH1F *child_pl_smear_pp = new TH1F("pl_smeared_pp",   "pL smeared pp",    500, 0, 500);
    TH1F *child_pl_smear_aa = new TH1F("pl_smeared_apap", "pL smeared apap",  500, 0, 500);
    TH1F *child_ps_smear_pp = new TH1F("ps_smeared_pp",   "pS+ smeared pp",   500, 0, 500);
    TH1F *child_ps_smear_aa = new TH1F("ps_smeared_apap", "pS+ smeared apap", 500, 0, 500);
    for (short nbin = 0; nbin < 500; ++nbin)
    {
	child_pl_smear_pp->SetBinContent(nbin + 1, decomp_pp->EvalSignalSmearedChild(0, child_pl_smear_pp->GetBinCenter(nbin + 1) + 1) + 1);
	child_pl_smear_aa->SetBinContent(nbin + 1, decomp_aa->EvalSignalSmearedChild(0, child_pl_smear_aa->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear_pp->SetBinContent(nbin + 1, decomp_pp->EvalSignalSmearedChild(1, child_ps_smear_pp->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear_aa->SetBinContent(nbin + 1, decomp_aa->EvalSignalSmearedChild(1, child_ps_smear_aa->GetBinCenter(nbin + 1) + 1) + 1);
    }

    TH1F **cfs_pp = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_pp[0] = static_cast<TH1F*>(input_cf_pp->Clone());
    cfs_pp[1] = child_pl_smear_pp;
    cfs_pp[2] = child_ps_smear_pp;

    TH1F **cfs_aa = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_aa[0] = static_cast<TH1F*>(input_cf_aa->Clone());
    cfs_aa[1] = child_pl_smear_aa;
    cfs_aa[2] = child_ps_smear_aa;

    TH1F *corrected_cf_pp = (TH1F*) input_cf_pp->Clone("CF_corrected_pp");
    corrected_cf_pp->Reset();
    get_corrected_cf(cfs_pp, lambda_pars_pp_th1, bsl_pars_pp, &corrected_cf_pp);

    TH1F *corrected_cf_aa = (TH1F*) input_cf_aa->Clone("CF_corrected_aa");
    corrected_cf_aa->Reset();
    get_corrected_cf(cfs_aa, lambda_pars_aa_th1, bsl_pars_aa, &corrected_cf_aa);

    TH1F **pot_histos = new TH1F*[5];
    std::vector<const char*> pot_names = {"1S0", "3P0", "3P1", "3P2", "1D2"};
    get_cats_potentials(&cats_pp, pot_histos, pot_names);

    printf("\n\e[1;36m  -->  \e[0;36mData Fit  \e[1;36m<--\e[0m\n\n");
    result.Print(std::cout); printf("\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Source Size  │   \e[0m%.2f ± %.2f\n", result.Parameter(10), result.Error(10));
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");

    TFile *ofile;
    create_output(ofile, "FitResults_", output_path, var);
    ofile->cd();

    input_cf_pp->Write("CF_pp");
    input_cf_aa->Write("CF_apap");

    corrected_cf_pp->Write("CF_corr_pp");
    corrected_cf_aa->Write("CF_corr_apap");

    input_me_orig_pp->Write("ME_orig_pp");
    input_me_orig_aa->Write("ME_orig_apap");

    fitter_pp->Write("fFitResult_pp");
    fitter_aa->Write("fFitResult_apap");

    baseline_pp->Write("fBl_pp");
    baseline_aa->Write("fBl_apap");

    genuine_pp->Write("fGenuine_pp");
    genuine_aa->Write("fGenuine_apap");

    genuine_smear_pp->Write("fGenuineSmeared_pp");
    genuine_smear_aa->Write("fGenuineSmeared_apap");

    child_pl_smear_pp->Write("pl_smeared_pp");
    child_pl_smear_aa->Write("pl_smeared_apap");
    child_ps_smear_pp->Write("ps_smeared_pp");
    child_ps_smear_aa->Write("ps_smeared_apap");

    pot_histos[0]->Write("hPotential_1S0");
    pot_histos[1]->Write("hPotential_3P0");
    pot_histos[2]->Write("hPotential_3P1");
    pot_histos[3]->Write("hPotential_3P2");
    pot_histos[4]->Write("hPotential_1D2");

    TDirectory *tdir_lam = ofile->mkdir("lam_pars");
    tdir_lam->cd();

    lambda_pars_pp_th1[0]->Write("genuine_pp");
    lambda_pars_pp_th1[1]->Write("lambda_pp");
    lambda_pars_pp_th1[2]->Write("sigma_pp");
    lambda_pars_pp_th1[3]->Write("flat_pp");
    lambda_pars_pp_th1[4]->Write("fake_pp");

    lambda_pars_aa_th1[0]->Write("genuine_apap");
    lambda_pars_aa_th1[1]->Write("lambda_apap");
    lambda_pars_aa_th1[2]->Write("sigma_apap");
    lambda_pars_aa_th1[3]->Write("flat_apap");
    lambda_pars_aa_th1[4]->Write("fake_apap");

    TDirectory *tdir_qa = ofile->mkdir("qa_plots");
    tdir_qa->cd();

    qa_plots->lam->Write();
    qa_plots->fmr->Write();
    qa_plots->fr->Write();
    qa_plots->bsl->Write();
    qa_plots->smear->Write();

    if (var->rsm)
    {
	qa_plots->frac->Write();
	qa_plots->mass->Write();
    }

    ofile->Close();
}

void cf_combined_fitter_pl(TString project_path, TString output_path, TString file_name_pp, TString file_name_aa, VAR *var)
{
    printf("\n  \e[1;36m-->  \e[0;36mEntering the combined fitter!  \e[1;36m<--\e[0m\n\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Input File   │\e[0m   %s\e[0m\n", file_name_pp.Data());
    printf("\e[1;34m  └───────────────┘\e[0m\n");

    VAR_RSM *var_rsm = new VAR_RSM(PROTON, LAMBDA, var->rsm);
    VAR_FMR *range_femto_pp = new VAR_FMR(var->system, var->fmr);
    VAR_FMR *range_femto_aa = new VAR_FMR(var->system, var->fmr);
    VAR_FR  *range_fit_pp = new VAR_FR(var->fr);
    VAR_FR  *range_fit_aa = new VAR_FR(var->fr);
    VAR_QA  *qa_plots = new VAR_QA(var);

    range_femto_pp->Max = 400;
    range_femto_pp->Min = 0;
    range_femto_aa->Max = 400;
    range_femto_aa->Min = 0;

    range_fit_pp->Max = 500;
    range_fit_pp->Min = 0;
    range_fit_aa->Max = 500;
    range_fit_aa->Min = 0;

    print_info_2(var, range_femto_pp, range_fit_pp);

    // fitter parameter entry list, i.e. the last parameter is shared between both fitters
    int fit_pars_pp[6] = { 0, 1, 2, 3, 4, 10 };	    // N, a, b, c, d, r
    int fit_pars_aa[6] = { 5, 6, 7, 8, 9, 10 };
    double global_pars[11] = {0.98, 0, 0, 0, 0, 0.98, 0, 0, 0, 0, 1};	    // init parameters

    std::unique_ptr<TH1F> input_cf_pp, input_me_pp, input_me_orig_pp;
    std::unique_ptr<TH1F> input_cf_aa, input_me_aa, input_me_orig_aa;
    get_input_pl(file_name_pp, input_cf_pp, input_me_pp, input_me_orig_pp, var, PP);
    get_input_pl(file_name_aa, input_cf_aa, input_me_aa, input_me_orig_aa, var, APAP);

    DLM_CommonAnaFunctions cats_setupper;
    cats_setupper.SetCatsFilesFolder(project_path + "input/filesCATS");
    TH2F *matrix_pl_res = cats_setupper.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
    TH2F *matrix_pl_dec = cats_setupper.GetResidualMatrix("pp", "pLambda");
    TH2F *matrix_ls0_dec = cats_setupper.GetResidualMatrix("pLambda", "pSigma0");
    TH2F *matrix_lxm_dec = cats_setupper.GetResidualMatrix("pLambda", "pXim");

    /* CATS Objects */
    CATS cats_pl_pp, cats_ps0_pp, cats_px0_pp, cats_pxm_pp;
    setup_cats(&cats_setupper, &cats_pl_pp,  range_femto_pp, "pl");
    setup_cats(&cats_setupper, &cats_ps0_pp, range_femto_pp, "ps0");
    setup_cats(&cats_setupper, &cats_px0_pp, range_femto_pp, "px0");
    setup_cats(&cats_setupper, &cats_pxm_pp, range_femto_pp, "pxm");

    CATS cats_pl_aa, cats_ps0_aa, cats_px0_aa, cats_pxm_aa;
    setup_cats(&cats_setupper, &cats_pl_aa,  range_femto_aa, "pl");
    setup_cats(&cats_setupper, &cats_ps0_aa, range_femto_aa, "ps0");
    setup_cats(&cats_setupper, &cats_px0_aa, range_femto_aa, "px0");
    setup_cats(&cats_setupper, &cats_pxm_aa, range_femto_aa, "pxm");

    /* Decomp Objects */
    DLM_Ck *ck_pl_pp, *ck_ps0_pp, *ck_px0_pp, *ck_pxm_pp;
    DLM_CkDecomposition *decomp_pl_pp, *decomp_ps0_pp, *decomp_px0_pp, *decomp_pxm_pp;
    setup_decomp(decomp_pl_pp,  ck_pl_pp,  &cats_pl_pp,  matrix_pl_res, range_femto_pp, "pl");
    setup_decomp(decomp_ps0_pp, ck_ps0_pp, &cats_ps0_pp, NULL,		range_femto_pp, "ps0");
    setup_decomp(decomp_px0_pp, ck_px0_pp, &cats_px0_pp, NULL,		range_femto_pp, "px0");
    setup_decomp(decomp_pxm_pp, ck_pxm_pp, &cats_pxm_pp, NULL,		range_femto_pp, "pxm");

    DLM_Ck *ck_pl_aa, *ck_ps0_aa, *ck_px0_aa, *ck_pxm_aa;
    DLM_CkDecomposition *decomp_pl_aa, *decomp_ps0_aa, *decomp_px0_aa, *decomp_pxm_aa;
    setup_decomp(decomp_pl_aa,  ck_pl_aa,  &cats_pl_aa,  matrix_pl_res, range_femto_aa, "pl");
    setup_decomp(decomp_ps0_aa, ck_ps0_aa, &cats_ps0_aa, NULL,		range_femto_aa, "ps0");
    setup_decomp(decomp_px0_aa, ck_px0_aa, &cats_px0_aa, NULL,		range_femto_aa, "px0");
    setup_decomp(decomp_pxm_aa, ck_pxm_aa, &cats_pxm_aa, NULL,		range_femto_aa, "pxm");

    decomp_pl_pp->AddContribution(0, LAM_PL["ps0"],  DLM_CkDecomp::cFeedDown, decomp_ps0_pp, matrix_ls0_dec);
    decomp_pl_pp->AddContribution(1, LAM_PL["xim"],  DLM_CkDecomp::cFeedDown, decomp_pxm_pp, matrix_lxm_dec);
    decomp_pl_pp->AddContribution(2, LAM_PL["flat"], DLM_CkDecomp::cFeedDown);
    decomp_pl_pp->AddContribution(3, LAM_PL["fake"], DLM_CkDecomp::cFake);

    decomp_pl_aa->AddContribution(0, LAM_PL["ps0"],  DLM_CkDecomp::cFeedDown, decomp_ps0_aa, matrix_ls0_dec);
    decomp_pl_aa->AddContribution(1, LAM_PL["xim"],  DLM_CkDecomp::cFeedDown, decomp_pxm_aa, matrix_lxm_dec);
    decomp_pl_aa->AddContribution(2, LAM_PL["flat"], DLM_CkDecomp::cFeedDown);
    decomp_pl_aa->AddContribution(3, LAM_PL["fake"], DLM_CkDecomp::cFake);

    // update CkDec objects
    decomp_pl_pp->Update(true, true);
    decomp_ps0_pp->Update(true, true);
    decomp_px0_pp->Update(true, true);
    decomp_pxm_pp->Update(true, true);

    decomp_pl_aa->Update(true, true);
    decomp_ps0_aa->Update(true, true);
    decomp_px0_aa->Update(true, true);
    decomp_pxm_aa->Update(true, true);

    // Fit function
    plFitFunctionSimple fitfun_pp(decomp_pl_pp);
    plFitFunctionSimple fitfun_aa(decomp_pl_aa);
    TF1 *fitter_pp = new TF1("Fit pp",	 fitfun_pp, range_fit_pp->Min, range_fit_pp->Max, 6);
    TF1 *fitter_aa = new TF1("Fit apap", fitfun_aa, range_fit_aa->Min, range_fit_aa->Max, 6);

    setup_fitter(fitter_pp, var, INIT_RADIUS[var->rsm][var->mult][var->mt]);
    setup_fitter(fitter_aa, var, INIT_RADIUS[var->rsm][var->mult][var->mt]);

    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange opt_range_pp, opt_range_aa;
    opt_range_pp.SetRange(range_fit_pp->Min, range_fit_pp->Max);
    opt_range_aa.SetRange(range_fit_aa->Min, range_fit_aa->Max);

    ROOT::Fit::BinData wrap_data_pp(opt, opt_range_pp);
    ROOT::Fit::BinData wrap_data_aa(opt, opt_range_aa);

    ROOT::Math::WrappedMultiTF1 wrap_fit_pp(*fitter_pp, 1);
    ROOT::Math::WrappedMultiTF1 wrap_fit_aa(*fitter_aa, 1);

    ROOT::Fit::FillData(wrap_data_pp, input_cf_pp.get());
    ROOT::Fit::FillData(wrap_data_aa, input_cf_aa.get());

    ROOT::Fit::Chi2Function chi2_pp(wrap_data_pp, wrap_fit_pp);
    ROOT::Fit::Chi2Function chi2_aa(wrap_data_aa, wrap_fit_aa);

    GlobalChi2 global_chi2(chi2_pp, chi2_aa, fit_pars_pp, fit_pars_aa);

    ROOT::Fit::Fitter fitter;
    setup_global_fitter(&fitter, var, global_pars);

    fitter.FitFCN(11, global_chi2, nullptr, wrap_data_pp.Size() + wrap_data_aa.Size(), true);
    ROOT::Fit::FitResult result = fitter.Result();

    fitter_pp->SetFitResult(result, fit_pars_pp);
    fitter_aa->SetFitResult(result, fit_pars_aa);

    double temp, total_chi2 = 0.;
    for (size_t nbin = input_cf_pp->FindBin(50); input_cf_pp->GetBinCenter(nbin) < 140; ++nbin)
    {
	temp = (std::pow(input_cf_pp->GetBinContent(nbin) - fitter_pp->Eval(input_cf_pp->GetBinCenter(nbin)), 2)) / std::pow(input_cf_pp->GetBinError(nbin), 2);
	total_chi2 += temp;
    }
    //std::cout << "pp fit Chi2 = " << total_chi2 << std::endl;
    out("pp fit Chi2 = " << total_chi2);

    fitter_pp->SetRange(opt_range_pp().first, opt_range_pp().second);
    fitter_aa->SetRange(opt_range_aa().first, opt_range_aa().second);

    double bsl_pars_pp[5], bsl_pars_aa[5];
    TF1 *baseline_pp = new TF1("Baseline pp",   bslfun, range_fit_pp->Min, range_fit_pp->Max, 5);
    TF1 *baseline_aa = new TF1("Baseline apap", bslfun, range_fit_aa->Min, range_fit_aa->Max, 5);
    for (size_t npar = 0; npar < 5; ++npar)
    {
	bsl_pars_pp[npar] = fitter_pp->GetParameter(npar);
	bsl_pars_aa[npar] = fitter_aa->GetParameter(npar);
	baseline_pp->SetParameter(npar, fitter_pp->GetParameter(npar));
	baseline_aa->SetParameter(npar, fitter_aa->GetParameter(npar));
    }

    GenuineFun genfun_pp(decomp_pl_pp), genfun_aa(decomp_pl_aa);
    TF1* genuine_pp = new TF1("Genuine pp", genfun_pp, range_fit_pp->Min, range_fit_pp->Max, 1);
    TF1* genuine_aa = new TF1("Genuine apap", genfun_aa, range_fit_aa->Min, range_fit_aa->Max, 1);
    genuine_aa->SetParameter(0, 1.);
    genuine_pp->SetParameter(0, 1.);

    GenuineFunSmear genfun_smear_pp(decomp_pl_pp), genfun_smear_aa(decomp_pl_aa);
    TF1* genuine_smear_pp = new TF1("Genuine Smeared pp", genfun_smear_pp, range_fit_pp->Min, range_fit_pp->Max, 1);
    TF1* genuine_smear_aa = new TF1("Genuine Smeared apap", genfun_smear_aa, range_fit_aa->Min, range_fit_aa->Max, 1);
    genuine_smear_pp->SetParameter(0, 1.);
    genuine_smear_aa->SetParameter(0, 1.);

    decomp_pl_pp->Update(true, true);
    decomp_pl_aa->Update(true, true);
    TH1F *child_pl_smear_pp = new TH1F("ps0_smeared_pp",   "pS0 smeared pp",    500, 0, 500);
    TH1F *child_pl_smear_aa = new TH1F("ps0_smeared_apap", "pS0 smeared apap",  500, 0, 500);
    TH1F *child_ps_smear_pp = new TH1F("pxm_smeared_pp",   "pS- smeared pp",   500, 0, 500);
    TH1F *child_ps_smear_aa = new TH1F("pxm_smeared_apap", "pS- smeared apap", 500, 0, 500);
    for (short nbin = 0; nbin < 500; ++nbin)
    {
	child_pl_smear_pp->SetBinContent(nbin + 1, decomp_pl_pp->EvalSignalSmearedChild(0, child_pl_smear_pp->GetBinCenter(nbin + 1) + 1) + 1);
	child_pl_smear_aa->SetBinContent(nbin + 1, decomp_pl_aa->EvalSignalSmearedChild(0, child_pl_smear_aa->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear_pp->SetBinContent(nbin + 1, decomp_pl_pp->EvalSignalSmearedChild(1, child_ps_smear_pp->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear_aa->SetBinContent(nbin + 1, decomp_pl_aa->EvalSignalSmearedChild(1, child_ps_smear_aa->GetBinCenter(nbin + 1) + 1) + 1);
    }

    TH1F **cfs_pp = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_pp[0] = static_cast<TH1F*>(input_cf_pp->Clone());
    cfs_pp[1] = child_pl_smear_pp;
    cfs_pp[2] = child_ps_smear_pp;

    TH1F **cfs_aa = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_aa[0] = static_cast<TH1F*>(input_cf_aa->Clone());
    cfs_aa[1] = child_pl_smear_aa;
    cfs_aa[2] = child_ps_smear_aa;

    printf("\n\e[1;36m  -->  \e[0;36mData Fit  \e[1;36m<--\e[0m\n\n");
    result.Print(std::cout); printf("\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Source Size  │   \e[0m%.2f ± %.2f\n", result.Parameter(10), result.Error(10));
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");

    TFile *ofile;
    create_output(ofile, "FitResults_", output_path, var);
    ofile->cd();

    input_cf_pp->Write("CF_pp");
    input_cf_aa->Write("CF_apap");

    input_me_orig_pp->Write("ME_orig_pp");
    input_me_orig_aa->Write("ME_orig_apap");

    fitter_pp->Write("fFitResult_pp");
    fitter_aa->Write("fFitResult_apap");

    baseline_pp->Write("fBl_pp");
    baseline_aa->Write("fBl_apap");

    genuine_pp->Write("fGenuine_pp");
    genuine_aa->Write("fGenuine_apap");

    genuine_smear_pp->Write("fGenuineSmeared_pp");
    genuine_smear_aa->Write("fGenuineSmeared_apap");

    child_pl_smear_pp->Write("ps0_smeared_pp");
    child_pl_smear_aa->Write("ps0_smeared_apap");
    child_ps_smear_pp->Write("psm_smeared_pp");
    child_ps_smear_aa->Write("psm_smeared_apap");

    TDirectory *tdir_qa = ofile->mkdir("qa_plots");
    tdir_qa->cd();

    qa_plots->lam->Write();
    qa_plots->fmr->Write();
    qa_plots->fr->Write();
    qa_plots->bsl->Write();
    qa_plots->smear->Write();

    if (var->rsm)
    {
	qa_plots->frac->Write();
	qa_plots->mass->Write();
    }

    ofile->Close();
}

