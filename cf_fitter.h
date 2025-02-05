#ifndef CF_FITTER_H
#define CF_FITTER_H

#include <map>

#include "TString.h"
#include "TH1F.h"
#include "TRandom3.h"

#define BINS(X) ((FEMTO_RANGE[X][1] - FEMTO_RANGE[X][0]) / 8.)

enum FITTER_TYPE { FIT_CATS = 0, FIT_RSM, FIT_BOTH };
enum CHARGE_TYPE { PP = 0, APAP, BOTH, PL, APAL, BOTHPL };

static double RSM_FRAC[3] = {0.6422, 0.90*RSM_FRAC[0], 1.10*RSM_FRAC[0]};
static double RSM_MASS[3] = {  1362, 0.90*RSM_MASS[0], 1.10*RSM_MASS[0]};
static double RSM_TAU[3]  = {  1.65, 0.90*RSM_TAU[0],  1.10*RSM_TAU[0]};

static const TString CHARGE_STR[] = {"pp", "apap", "combined", "pl", "apal", "combined"};

static double
FEMTO_RANGE[][3] = {
    {0, 280, BINS(0)},
    {0, 240, BINS(1)},
    {0, 320, BINS(2)},
};

static double
FEMTO_RANGE_PL[][3] = {
    {0, 400, BINS(0)},
    {0, 360, BINS(1)},
    {0, 440, BINS(2)},
};

static double
FIT_RANGE[][3] = {
    {0, 400},
    {0, 360},
    {0, 440},
};

static std::map<std::string, double>
LAM_PL = {
    {"gen",  0.4969},
    {"ps0",  0.1677},
    {"xi0",  0.0000},
    {"xim",  0.0832},
    {"flat", 0.2038},
    {"fake", 0.0422},
};

static double
INIT_RADIUS[2][3][8] = {
    {
	{0, 1.30, 1.25, 1.15, 1.10, 1.05, 0.95, 0.80},
	{0, 1.20, 1.15, 1.10, 1.05, 1.00, 0.90, 0.70},
	{0, 1.10, 1.05, 1.00, 0.95, 0.90, 0.85, 0.60}
    },{
	{0, 1.10, 1.05, 1.00, 0.95, 0.90, 0.80, 0.70},
	{0, 1.00, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60},
	{0, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.50}
    }
};

typedef struct
VAR {
    int system = 0;
    int charge = 0;
    int sample = 0;
    int stat   = 0;
    int rsm    = 0;
    int mt     = 2;
    int mult   = 0;
    int lam    = 0;
    int fmr    = 0;
    int fr     = 0;
    int bsl    = 0;
    int smear  = 0;
    int frac   = 0;
    int mass   = 0;
    int normal = 0;
} VAR;

inline void
set_vars(int argc, char *argv[], VAR *var)
{
    if (argc > 3) var->charge = atoi(argv[3]);
    if (argc > 4) var->sample = atoi(argv[4]);
    if (argc > 5) var->stat   = atoi(argv[5]);
    if (argc > 6) var->rsm    = atoi(argv[6]);
    if (argc > 7) var->mt     = atoi(argv[7]);
    if (argc > 8) var->mult   = atoi(argv[8]);

    if (!var->stat)
    {
	TRandom3 *random = new TRandom3(var->sample);
	if (var->sample)
	{
	    var->lam   = random->Integer(3);
	    var->fmr   = random->Integer(3);
	    var->fr    = random->Integer(3);
	    //var->bsl   = random->Integer(2);
	    var->bsl   = 0;
	    //var->smear = random->Integer(2);
	    var->smear = 0;
	    var->frac  = random->Integer(3);
	    var->mass  = random->Integer(3);
	}

	if (argc >  9) var->lam   = atoi(argv[9]);
	if (argc > 10) var->fmr   = atoi(argv[10]);
	if (argc > 11) var->fr    = atoi(argv[11]);
	if (argc > 12) var->bsl   = atoi(argv[12]);
	if (argc > 13) var->smear = atoi(argv[13]);
	if (argc > 14) var->frac  = atoi(argv[14]);
	if (argc > 15) var->mass  = atoi(argv[15]);
    }
}

typedef struct
VAR_QA {
    TH1F *lam   = new TH1F("qa_lambda",	"qa_lambda", 10, -0.5, 9.5);
    TH1F *fmr   = new TH1F("qa_femto",	"qa_femto",  10, -0.5, 9.5);
    TH1F *fr    = new TH1F("qa_fit",	"qa_fit",    10, -0.5, 9.5);
    TH1F *bsl   = new TH1F("qa_bsl",	"qa_bsl",    10, -0.5, 9.5);
    TH1F *smear = new TH1F("qa_smear",	"qa_smear",  10, -0.5, 9.5);
    TH1F *frac  = new TH1F("qa_frac",	"qa_frac",   10, -0.5, 9.5);
    TH1F *mass  = new TH1F("qa_mass",	"qa_mass",   10, -0.5, 9.5);

    VAR_QA(VAR *var)
    {
	lam->Fill(var->lam);
	fmr->Fill(var->fmr);
	fr->Fill(var->fr);
	bsl->Fill(var->bsl);
	smear->Fill(var->smear);
	frac->Fill(var->frac);
	mass->Fill(var->mass);
    }

    ~VAR_QA()
    {
	delete lam;
	delete fmr;
	delete fr;
	delete bsl;
	delete smear;
	delete frac;
	delete mass;
    }
} VAR_QA;

typedef struct
VAR_FMR {
    double Min, Max, Bins;
    VAR_FMR() : Min(FEMTO_RANGE[0][0]), Max(FEMTO_RANGE[0][1]), Bins(FEMTO_RANGE[0][2]) {}
    VAR_FMR(int fmr) :
	Min(FEMTO_RANGE[fmr][0]),
	Max(FEMTO_RANGE[fmr][1]),
	Bins(FEMTO_RANGE[fmr][2])
    {}
} VAR_FMR;

typedef struct
VAR_FMR_PL {
    double Min, Max, Bins;
    VAR_FMR_PL() : Min(FEMTO_RANGE_PL[0][0]), Max(FEMTO_RANGE_PL[0][1]), Bins(FEMTO_RANGE_PL[0][2]) {}
    VAR_FMR_PL(int fmr) :
	Min(FEMTO_RANGE_PL[fmr][0]),
	Max(FEMTO_RANGE_PL[fmr][1]),
	Bins(FEMTO_RANGE_PL[fmr][2])
    {}
} VAR_FMR_PL;

typedef struct
VAR_FR {
    double Min, Max;
    VAR_FR() : Min(FIT_RANGE[0][0]), Max(FIT_RANGE[0][1]) {}
    VAR_FR(int fr) : Min(FIT_RANGE[fr][0]), Max(FIT_RANGE[fr][1]) {}
} VAR_FR;

typedef struct
VAR_RSM {
    double frac, mass, tau;
    VAR_RSM() : frac(RSM_FRAC[0]), mass(RSM_MASS[0]), tau(RSM_TAU[0]) {}
    VAR_RSM(int rsm) : frac(RSM_FRAC[rsm]), mass(RSM_MASS[rsm]), tau(RSM_TAU[0]) {}
} VAR_RSM;

void cf_fitter(TString, TString, TString, VAR*);
void cf_fitter_pl(TString, TString, TString, VAR*);
void cf_combined_fitter(TString, TString, TString, TString, TString, VAR*);
void cf_combined_fitter_pl(TString, TString, TString, TString, VAR*);

#endif
