/***********************************************************************
/
/  COSMOLOGY: COMPUTES THE TIME (IN CODE UNITS) FROM THE GIVEN REDSHIFT
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"
 
// function prototypes

int CosmologyTableComputeTimeFromRedshift(FLOAT z, FLOAT *time);


int CosmologyComputeTimeFromRedshift(FLOAT Redshift, FLOAT *TimeCodeUnits)
{

  float TimeUnits = 2.52e17/sqrt(OmegaMatterNow)/HubbleConstantNow/
    POW(1 + InitialRedshift, FLOAT(1.5));

  if (UseModifiedLambda == TRUE) {

    FLOAT a = 1 / (1 + Redshift);

    if (ABS(VacuumAlpha) < 1e-10) {
      ENZO_FAIL("VacuumAlpha cannot be zero.\n");
    }

    //FLOAT tau = TimeHubble0 * sqrt(OmegaLambdaNow) / ABS(VacuumAlpha);
    FLOAT aref = POW((OmegaMatterNow / (ABS(VacuumAlpha) * OmegaLambdaNow)), 1./3.);
    FLOAT tH = 1 / (HubbleConstantNow*3.24e-18);
    FLOAT tref = POW((ABS(VacuumAlpha) * OmegaLambdaNow),-0.5) * tH;
    FLOAT bigA = a / aref;
    FLOAT tau;

    if (VacuumAlpha < 0) {
      tau = (2./3.) * asin(POW(bigA, 1.5));
      if (BeforeTurnAround == TRUE) {
        tau = 2. * pi / 3. - tau;
      }
    }
    else {
      tau = (2./3.) * asinh(POW(bigA, 1.5));
    }

    *TimeCodeUnits = tau * tref / TimeUnits;
  }

  else {
 
  FLOAT TimeHubble0;
 
  /* Find Omega due to curvature. */
 
  float OmegaCurvatureNow = 1 - OmegaMatterNow -
    OmegaLambdaNow - OmegaRadiationNow;

  /* Interpolate from a(t) table. */
 
  if (CosmologyTableComputeTimeFromRedshift(Redshift, &TimeHubble0) == FAIL) {
    ENZO_FAIL("Error in CosmologyTableComputeTime.\n");
  }

  /* Now convert from Time * H0 to code units (see also CosmologyGetUnits). */

  *TimeCodeUnits = TimeHubble0 / (HubbleConstantNow*3.24e-18) / TimeUnits;

  }
 
  return SUCCESS;
}
