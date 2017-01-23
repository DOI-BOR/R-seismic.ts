#include "resp.h"
#include "resp_proto.h"

/* be sure to include these last so they can override some windows defs */
#include <R.h>
#include <Rinternals.h>

/* resp() wrapper for R */
DllExport SEXP CALLresp(SEXP in_type_s, SEXP out_type_s,
	SEXP ts_d, SEXP dt_d, SEXP ptap_d,
	SEXP lambda_d, SEXP tau_range_d, SEXP tau_d,
	SEXP rs_method_i, SEXP tau_si_range_d)
{
	int ii, jj, pcnt = 0, ndamp = -1, have_damping = FALSE, nper = -1, have_periods = FALSE;
	double *damp = NULL, *per = NULL, *si = NULL, **rspect = NULL, *dp;
	char rstype[128], units[128], siunits[128];
	SEXP damping_d, periods_d, rspectra_d, rm_s, si_d, tausi_d, rs_type_s, units_s, si_units_s;	/* return values */
	SEXP ret_l, names_s, dimnames_s;

	char in_type = asChar(in_type_s) == NA_STRING ? 'a' : tolower(*CHAR(STRING_ELT(in_type_s,0)));
	char out_type = asChar(out_type_s) == NA_STRING ? 'a' : tolower(*CHAR(STRING_ELT(out_type_s,0)));
	double *ts = REAL(ts_d);
	int len = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	double ptap = ! R_FINITE(asReal(ptap_d)) ? -1 : REAL(ptap_d)[0];
	RespMethod rm = asInteger(rs_method_i) == NA_INTEGER ? RM_IIR : INTEGER(rs_method_i)[0];
	double tau_lo = ! R_FINITE(asReal(tau_range_d)) || length(tau_range_d) != 2 ? -1 : REAL(tau_range_d)[0];
	double tau_hi = ! R_FINITE(asReal(tau_range_d)) || length(tau_range_d) != 2 ? -1 : REAL(tau_range_d)[1];
	double tausi_lo = ! R_FINITE(asReal(tau_si_range_d)) || length(tau_si_range_d) != 2 ? -1 : REAL(tau_si_range_d)[0];
	double tausi_hi = ! R_FINITE(asReal(tau_si_range_d)) || length(tau_si_range_d) != 2 ? -1 : REAL(tau_si_range_d)[1];

	/* check for non-default input */
	if ( R_FINITE(asReal(lambda_d)) && length(lambda_d) > 0 ) {
		damp = REAL(lambda_d);
		ndamp = length(lambda_d);
		have_damping = TRUE;
	}
	if ( R_FINITE(asReal(tau_d)) && length(tau_d) > 0 ) {
		per = REAL(tau_d);
		nper = length(tau_d);
		have_periods = TRUE;
	}

	/* call resp */
	resp(in_type, out_type, ts, len, dt, ptap, &damp, &ndamp,
			tau_lo, tau_hi, &per, &nper,
			rm, &rspect, rstype, units,
			&si, &tausi_lo, &tausi_hi, siunits);

	/* allocate space for R structures for damping, periods, rspectra, and si */
	damping_d = PROTECT(allocVector(REALSXP, ndamp)); pcnt++;
	dp = REAL(damping_d);
	for ( ii = 0 ; ii < ndamp ; ii++ )
		dp[ii] = damp[ii];

	periods_d = PROTECT(allocVector(REALSXP, nper)); pcnt++;
	dp = REAL(periods_d);
	for ( ii = 0 ; ii < nper ; ii++ )
		dp[ii] = per[ii];

	rspectra_d = PROTECT(allocMatrix(REALSXP, nper, ndamp)); pcnt++;
	dp = REAL(rspectra_d);	/* 1-D pointer into matrix */
	/* copy response spectra values to R data structure */
	for ( ii = 0 ; ii < ndamp ; ii++ )
		for ( jj = 0 ; jj < nper ; jj++ )
			dp[ii * nper + jj] = (rspect[ii])[jj];

	/* set the dim names */
	dimnames_s = PROTECT(allocVector(VECSXP, 2)); pcnt++;
	names_s = PROTECT(allocVector(VECSXP, nper)); pcnt++;
	for ( ii = 0 ; ii < nper ; ii++ ) {
		char label[32];
		sprintf(label,"%.3f sec.",per[ii]);
		SET_VECTOR_ELT(names_s, ii, mkChar(label));
	}
	SET_VECTOR_ELT(dimnames_s, 0, names_s);
	names_s = PROTECT(allocVector(VECSXP, ndamp)); pcnt++;
	for ( ii = 0 ; ii < ndamp ; ii++ ) {
		char label[32];
		sprintf(label,"%.1f pct.",damp[ii]*100);
		SET_VECTOR_ELT(names_s, ii, mkChar(label));
	}
	SET_VECTOR_ELT(dimnames_s, 1, names_s);
	setAttrib(rspectra_d, R_DimNamesSymbol, dimnames_s);

	si_d = PROTECT(allocVector(REALSXP, ndamp)); pcnt++;
	dp = REAL(si_d);
	for ( ii = 0 ; ii < ndamp ; ii++ )
		dp[ii] = si[ii];
	if ( ndamp > 1 )
		setAttrib(si_d, R_NamesSymbol, names_s);

	tausi_d = PROTECT(allocVector(REALSXP, 2)); pcnt++;
	REAL(tausi_d)[0] = tausi_lo;
	REAL(tausi_d)[1] = tausi_hi;

	rs_type_s = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(rs_type_s, 0, mkChar(rstype));
	units_s = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(units_s, 0, mkChar(units));
	si_units_s = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(si_units_s, 0, mkChar(siunits));
	rm_s = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(rm_s, 0, mkChar(rm == RM_IIR ? "IIR" :
										rm == RM_BZ ? "Bi-linear Z-transform" :
										rm == RM_LS ? "Niagam and Jennings (1969)" :
										rm == RM_CONV ? "Circular convolution" : "Unknown"));

	/* put the return values into a list */
	ret_l = PROTECT(allocVector(VECSXP, 9)); pcnt++;
	names_s = PROTECT(allocVector(VECSXP, 9)); pcnt++;
	SET_VECTOR_ELT(ret_l, 0, damping_d);
	SET_VECTOR_ELT(names_s, 0, mkChar("damping"));
	SET_VECTOR_ELT(ret_l, 1, periods_d);
	SET_VECTOR_ELT(names_s, 1, mkChar("periods"));
	SET_VECTOR_ELT(ret_l, 2, rspectra_d);
	SET_VECTOR_ELT(names_s, 2, mkChar("rspect"));
	SET_VECTOR_ELT(ret_l, 3, rs_type_s);
	SET_VECTOR_ELT(names_s, 3, mkChar("rs.type"));
	SET_VECTOR_ELT(ret_l, 4, units_s);
	SET_VECTOR_ELT(names_s, 4, mkChar("rs.units"));
	SET_VECTOR_ELT(ret_l, 5, rm_s);
	SET_VECTOR_ELT(names_s, 5, mkChar("rs.method"));
	SET_VECTOR_ELT(ret_l, 6, si_d);
	SET_VECTOR_ELT(names_s, 6, mkChar("SI"));
	SET_VECTOR_ELT(ret_l, 7, tausi_d);
	SET_VECTOR_ELT(names_s, 7, mkChar("SI.per.range"));
	SET_VECTOR_ELT(ret_l, 8, si_units_s);
	SET_VECTOR_ELT(names_s, 8, mkChar("SI.units"));
	setAttrib(ret_l, R_NamesSymbol, names_s);

	/* free space created by resp */
	if ( ! have_damping )
		free(damp);
	if ( ! have_periods )
		free(per);
	for ( ii = 0 ; ii < ndamp ; ii ++ )
		if ( rspect[ii] != NULL )
			free(rspect[ii]);	/* free the response spectra arrays */
	free(rspect);	/* free the pointer array to the response spectra */
	if ( ndamp > 0 )
		free(si);

	UNPROTECT(pcnt);

	return ret_l;
}
