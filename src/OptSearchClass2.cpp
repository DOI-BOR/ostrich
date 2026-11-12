/**************************************************************************************************************************************************************

One-dimensional line search used to compute the optimum step size along a given search direction. The minimum is first bracketed and then refined using
either the Golden Section or Brent method. Parameter sets are evaluated through a user-supplied evaluator so that the search can be driven by the asynchronous
MPI work-farming functionality provided by the Algorithm base class. Independent evaluations are submitted to the evaluator together so that they can be
solved concurrently.

**************************************************************************************************************************************************************/

// Include the class header
#include "OptSearchClass2.h"

// Golden section constants
static const double gCGOLD2 = (0.5 * (3.00 - sqrt(5.00)));
static const double gGOLD2 = (2.0 - (0.5 * (3.00 - sqrt(5.00))));

/**************************************************************************************************************************************************************
CTOR

Initializes the line search. All settings are assigned default values and then the user input file is checked for overriding values.
**************************************************************************************************************************************************************/
OptSearchClass2::OptSearchClass2(ParameterGroup* pParamGroup, Evaluator evaluator) {

    m_pParamGroup = pParamGroup;
    m_Evaluate = evaluator;
    m_NumParams = m_pParamGroup->GetNumParams();
    m_Dir = std::vector<double>(m_NumParams, 0.00);

    InitFromFile(GetInFileName());

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the line search.
**************************************************************************************************************************************************************/
void OptSearchClass2::Destroy(void) {

    m_StartPoint.clear();
    m_Dir.clear();

    IncDtorCount();
}

/**************************************************************************************************************************************************************
WriteMetrics()

Reports on the setup of the line search and also various run-time metrics.
**************************************************************************************************************************************************************/
void OptSearchClass2::WriteMetrics(FILE* pFile) {

    fprintf(pFile, "\nOne Dimensional Search Metrics\n");
    fprintf(pFile, "Search Type          : ");
    if (m_SearchType == SEARCH_GOLDEN_SECTION) { fprintf(pFile, "Golden Section\n"); }
    else { fprintf(pFile, "Brent\n"); }
    fprintf(pFile, "Convergence Val      : %E\n", m_SearchConvVal);
    fprintf(pFile, "Bound Min Evals      : %d\n", m_BoundMinCount);
    fprintf(pFile, "Golden Section Evals : %d\n", m_GoldSectCount);
    fprintf(pFile, "Brent Evals          : %d\n", m_BrentCount);
}

/**************************************************************************************************************************************************************
InitFromFile()

Reads the line search configuration from the input file.
**************************************************************************************************************************************************************/
void OptSearchClass2::InitFromFile(IroncladString pFileName) {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];

    // Read in the configuration
    pFile = fopen(pFileName, "r");
    if (pFile == NULL) {
        // Couldn't open file, use defaults and log the error
        LogError(ERR_FILE_IO, "Couldn't open 1-D search config. file. Using Defaults");
        return;
    }

    // Make sure the correct tokens are present
    if (CheckToken(pFile, "Begin1dSearch", pFileName) == true) {
        // Identify the end of the block
        FindToken(pFile, "End1dSearch", pFileName);
        rewind(pFile);

        // Import data starting from the beginning of the block
        FindToken(pFile, "Begin1dSearch", pFileName);
        line = GetNxtDataLine(pFile, pFileName);

        while (strstr(line, "End1dSearch") == NULL) {
            if (strstr(line, "1dSearchConvergeVal") != NULL) {
                // Set the convergence value
                sscanf(line, "%s %lf", tmp, &m_SearchConvVal);
            }
            else if (strstr(line, "1dSearchMethod") != NULL) {
                // Set the search method
                MyStrLwr(line);
                if (strstr(line, "brent") != NULL) { m_SearchType = SEARCH_BRENT; }
                else if (strstr(line, "goldensection") != NULL) { m_SearchType = SEARCH_GOLDEN_SECTION; }
                else {
                    strcpy(tmp, "Unknown Search Method");
                    LogError(ERR_FILE_IO, tmp);
                }
            }
            else {
                sprintf(tmp, "Unknown token: %s", line);
                LogError(ERR_FILE_IO, tmp);
            }

            // Get the next line
            line = GetNxtDataLine(pFile, pFileName);
        }
    }

    fclose(pFile);
}

/**************************************************************************************************************************************************************
CreateAlphaPoint()

Creates the parameter set located a step size of <alpha> along the search direction from the start point. Any parameter that steps out of its bounds is
moved half the distance to the violated bound.
**************************************************************************************************************************************************************/
std::vector<double> OptSearchClass2::CreateAlphaPoint(double alpha) {

    std::vector<double> point = m_StartPoint;

    for (int entryParameter = 0; entryParameter < m_NumParams; entryParameter++) {
        double oldValue = point[entryParameter];
        point[entryParameter] = oldValue + alpha * m_Dir[entryParameter];

        // If it's out of bounds, move half the distance to upr/lwr
        double upr = m_pParamGroup->GetParamPtr(entryParameter)->GetUpperBoundTransformed();
        double lwr = m_pParamGroup->GetParamPtr(entryParameter)->GetLowerBoundTransformed();
        if (point[entryParameter] > upr) { point[entryParameter] = (upr + oldValue) / 2.00; }
        if (point[entryParameter] < lwr) { point[entryParameter] = (oldValue + lwr) / 2.00; }
    }

    return point;
}

/**************************************************************************************************************************************************************
CalcF()

Computes the objective function for each of the step sizes in <alphas>. The step sizes are evaluated together. If a better minimum than fmin is found, fmin
and xmin are replaced with the new minimum value and corresponding parameter values (in the order of the step sizes).
**************************************************************************************************************************************************************/
std::vector<double> OptSearchClass2::CalcF(std::vector<double> alphas, double& fmin, std::vector<double>& xmin) {

    // Create the parameter sets for each step size
    std::vector<std::vector<double>> points;
    for (int entryAlpha = 0; entryAlpha < (int)alphas.size(); entryAlpha++) {
        points.push_back(CreateAlphaPoint(alphas[entryAlpha]));
    }

    // Evaluate the parameter sets
    std::vector<double> objectives = m_Evaluate(points);

    // Update the optimal, if appropriate
    for (int entryAlpha = 0; entryAlpha < (int)alphas.size(); entryAlpha++) {
        if (objectives[entryAlpha] < fmin) {
            fmin = objectives[entryAlpha];
            xmin = points[entryAlpha];
        }
    }

    return objectives;
}

/**************************************************************************************************************************************************************
CalcF()

Computes the objective function for a single step size of <alpha>.
**************************************************************************************************************************************************************/
double OptSearchClass2::CalcF(double alpha, double& fmin, std::vector<double>& xmin) {

    std::vector<double> alphas = { alpha };
    return CalcF(alphas, fmin, xmin)[0];
}

/**************************************************************************************************************************************************************
CalcStepSize()

Computes the optimum step size in the direction given by <direction>, starting from <startPoint> (which has an objective of <startObjective>). If a better
minimum than fmin is found, fmin and xmin are replaced with the new minimum value and corresponding parameter values.
**************************************************************************************************************************************************************/
double OptSearchClass2::CalcStepSize(std::vector<double> startPoint, double startObjective, std::vector<double> direction, double& fmin,
                                     std::vector<double>& xmin) {

    // Store the current setting and direction
    m_StartPoint = startPoint;
    m_StartObjective = startObjective;
    m_Dir = direction;

    // Bound the minimum
    Bracket bracket = BracketMinimum(-1.00, 1.00, fmin, xmin);

    // Refine the minimum
    MinimumPoint minimum;
    if (m_SearchType == SEARCH_GOLDEN_SECTION) { minimum = GoldSect(bracket, fmin, xmin); }
    else { minimum = Brent(bracket, fmin, xmin); }

    return minimum.x;
}

/**************************************************************************************************************************************************************
BracketMinimum()

Brackets the minimum using a standard Golden Section algorithm. The two initial step sizes are evaluated together.

Returns a bracket containing points a, b and c along with f(a), f(b), and f(c). If a better minimum than fmin is found, fmin and xmin are replaced with the
new minimum value and corresponding parameter values.
**************************************************************************************************************************************************************/
OptSearchClass2::Bracket OptSearchClass2::BracketMinimum(double a, double b, double& fmin, std::vector<double>& xmin) {

    double ulim, u, r, q, fu, dum, denom;
    double fa, fb, c, fc;

    // Initial values
    double Fcur = m_StartObjective;

    // Evaluate the two initial step sizes together
    a = LimitStepSize(a);
    b = LimitStepSize(b);
    std::vector<double> initial = CalcF(std::vector<double>{ a, b }, fmin, xmin);
    fa = initial[0];
    fb = initial[1];
    m_BoundMinCount += 2;

    if (fb > fa) {
        // Switch roles of a and b
        dum = a; a = b; b = dum;
        dum = fa; fa = fb; fb = dum;
    }

    // Check to see if initial fb and fa already bracket a minimum
    if ((fb >= Fcur) && (fa >= Fcur)) {
        c = b;
        fc = fb;
        b = 0.00;
        fb = Fcur;
    }
    else {
        c = b + (gGOLD2 * (b - a));
        c = LimitStepSize(c);
        fc = CalcF(c, fmin, xmin);
        m_BoundMinCount++;

        while (fb > fc) {
            r = (b - a) * (fb - fc);
            q = (b - c) * (fb - fa);
            denom = MyMax(fabs(q - r), NEARLY_ZERO);
            if ((q - r) < 0) denom = -denom;
            u = (b) - ((b - c) * q - (b - a) * r) / (2.0 * denom);
            ulim = b + 100 * (c - b);

            // Test various possibilities
            if ((b - u) * (u - c) > 0.0) {
                u = LimitStepSize(u);
                fu = CalcF(u, fmin, xmin);
                m_BoundMinCount++;

                if (fu <= fc) {
                    // Minimum between b and c
                    a = b;
                    b = u;
                    fa = fb;
                    fb = fu;
                    break;
                }
                else if (fu >= fb) {
                    // Minimum between a and u
                    c = u;
                    fc = fu;
                    break;
                }
                u = (c) + gGOLD2 * (c - b);
                u = LimitStepSize(u);
                fu = CalcF(u, fmin, xmin);
                m_BoundMinCount++;
            }
            else if ((c - u) * (u - ulim) > 0.0) {
                u = LimitStepSize(u);
                fu = CalcF(u, fmin, xmin);
                m_BoundMinCount++;
                if (fu < fc) {
                    dum = c + gGOLD2 * (c - b);
                    b = c;
                    c = u;
                    u = dum;

                    dum = CalcF(u, fmin, xmin);
                    fb = fc;
                    fc = fu;
                    fu = dum;

                    m_BoundMinCount++;
                }
            }
            else if ((u - ulim) * (ulim - c) >= 0.0) {
                u = ulim;
                u = LimitStepSize(u);
                fu = CalcF(u, fmin, xmin);
                m_BoundMinCount++;
            }
            else {
                u = (c) + gGOLD2 * (c - b);
                u = LimitStepSize(u);
                fu = CalcF(u, fmin, xmin);
                m_BoundMinCount++;
            }
            a = b;
            b = c;
            c = u;

            fa = fb;
            fb = fc;
            fc = fu;
        }
    }

    Bracket bracket;
    bracket.a = a;
    bracket.b = b;
    bracket.c = c;
    bracket.fa = fa;
    bracket.fb = fb;
    bracket.fc = fc;

    return bracket;
}

/**************************************************************************************************************************************************************
GoldSect()

Calculates the minimum using the golden section method, as described by Vanderplaats (page 57-60). The two initial interior points are evaluated together.

Returns the location of the minimum point. If a better minimum than fmin is found, fmin and xmin are replaced with the new minimum value and corresponding
parameter values.
**************************************************************************************************************************************************************/
OptSearchClass2::MinimumPoint OptSearchClass2::GoldSect(Bracket bracket, double& fmin, std::vector<double>& xmin) {

    double x1, x2, F1, F2, xL, xU, Fbest, xbest;

    // Determine the initial minimum
    Fbest = bracket.fa; xbest = bracket.a;
    if (Fbest > bracket.fb) { Fbest = bracket.fb; xbest = bracket.b; }
    if (Fbest > bracket.fc) { Fbest = bracket.fc; xbest = bracket.c; }

    xL = bracket.a;
    xU = bracket.c;

    // Calculate the number of iterations required to reduce the initial width to the desired convergence value
    int its = 0;
    double width = fabs(xL - xU);
    if (width > m_SearchConvVal) {
        its = (int)((log10(m_SearchConvVal / width) / log10(1.00 - gCGOLD2)) + 1.00);
    }

    Write1dSearch2(WRITE_GSECT, its + 2);

    // Evaluate the two initial interior points together
    x1 = ((1.00 - gCGOLD2) * xL) + (gCGOLD2 * xU);
    x1 = LimitStepSize(x1);
    x2 = ((1.00 - gCGOLD2) * xU) + (gCGOLD2 * xL);
    x2 = LimitStepSize(x2);
    std::vector<double> initial = CalcF(std::vector<double>{ x1, x2 }, fmin, xmin);
    F1 = initial[0];
    F2 = initial[1];
    if (Fbest > F1) { Fbest = F1; xbest = x1; }
    if (Fbest > F2) { Fbest = F2; xbest = x2; }
    m_GoldSectCount += 2;
    Write1dSearch2(1, 0);
    Write1dSearch2(2, 0);

    // Golden Section iteration
    for (int entryIteration = 0; entryIteration < its; entryIteration++) {
        Write1dSearch2(entryIteration + 3, 0);
        if (F1 > F2) {
            xL = x1;
            x1 = x2;
            F1 = F2;
            x2 = ((1.00 - gCGOLD2) * xU) + (gCGOLD2 * xL);
            x2 = LimitStepSize(x2);
            F2 = CalcF(x2, fmin, xmin);
            if (Fbest > F2) { Fbest = F2; xbest = x2; }
            m_GoldSectCount++;
        }
        else {
            xU = x2;
            x2 = x1;
            F2 = F1;
            x1 = ((1.00 - gCGOLD2) * xL) + (gCGOLD2 * xU);
            x1 = LimitStepSize(x1);
            F1 = CalcF(x1, fmin, xmin);
            if (Fbest > F1) { Fbest = F1; xbest = x1; }
            m_GoldSectCount++;
        }
    }

    MinimumPoint minimum;
    minimum.x = xbest;
    minimum.fx = Fbest;

    Write1dSearch2(WRITE_ENDED, 0);
    return minimum;
}

/**************************************************************************************************************************************************************
Brent()

Finds the minimum, given the bracketing triplet, using the well-known Brent's method. Reverts to the Golden Section method if the maximum number of iterations
is exceeded or the search stalls on the first iteration.

Returns the minimum point and the value of the objective function at that point. If a better minimum than fmin is found, fmin and xmin are replaced with the
new minimum value and corresponding parameter values.
**************************************************************************************************************************************************************/
OptSearchClass2::MinimumPoint OptSearchClass2::Brent(Bracket bracket, double& fmin, std::vector<double>& xmin) {

    int iter;
    double ax, bx, cx, tol, xbest = 0.00, fbest = 0.00;
    double a, b, d = 0.00, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double e = 0.0;

    ax = bracket.a;
    bx = bracket.b;
    cx = bracket.c;
    tol = m_SearchConvVal;

    if (ax < cx) { a = ax; } else { a = cx; }
    if (ax > cx) { b = ax; } else { b = cx; }
    x = w = v = bx;
    fw = fv = fx = bracket.fb;

    Write1dSearch2(WRITE_BRENT, m_MaxBrentIters);
    for (iter = 1; iter <= m_MaxBrentIters; iter++) {
        Write1dSearch2(iter, 0);
        xbest = x;
        fbest = fx;

        xm = 0.5 * (a + b);
        tol1 = tol * fabs(x) + NEARLY_ZERO;
        tol2 = 2.0 * tol1;
        if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
            // Converged on minimum, so break out of the loop
            break;
        }
        if (fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) { p = -p; }
            q = fabs(q);
            etemp = e;
            e = d;
            if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                if (x >= xm) { e = (a - x); } else { e = (b - x); }
                d = gCGOLD2 * e;
            }
            else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2) {
                    if ((xm - x) < 0) d = -fabs(tol1);
                    else d = fabs(tol1);
                }
            }
        }
        else {
            if (x >= xm) { e = (a - x); } else { e = (b - x); }
            d = gCGOLD2 * e;
        }
        if (fabs(d) >= tol1) {
            u = x + d;
        }
        else {
            if (d >= 0) { u = x + tol1; }
            else { u = x - tol1; }
        }

        // Prevent step size from exceeding parameter bounds
        u = LimitStepSize(u);
        fu = CalcF(u, fmin, xmin);
        m_BrentCount++;

        if (fu == fx) { break; } // Multi-modal?
        if (fu < fx) {
            // Found a better minimum
            if (u >= x) { a = x; } else { b = x; }
            v = w;
            w = x;
            x = u;

            fv = fw;
            fw = fx;
            fx = fu;
        }
        else {
            if (u < x) { a = u; } else { b = u; }
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }

    // If we exceed the max. iterations or stalled at first iteration, revert to the trusty Golden Section method
    if ((iter > m_MaxBrentIters) || (iter == 1)) {
        Write1dSearch2(WRITE_SWTCH, 0);
        return GoldSect(bracket, fmin, xmin);
    }

    Write1dSearch2(WRITE_ENDED, 0);
    MinimumPoint minimum;
    minimum.x = xbest;
    minimum.fx = fbest;
    return minimum;
}

/**************************************************************************************************************************************************************
LimitStepSize()

Prevents the step size from exceeding parameter bounds.

Returns the largest allowable step size given the parameter bounds and the search direction.
**************************************************************************************************************************************************************/
double OptSearchClass2::LimitStepSize(double alpha) {

    // Check each design parameter to see if it can accommodate the given step size
    for (int entryParameter = 0; entryParameter < m_NumParams; entryParameter++) {
        double opi = m_StartPoint[entryParameter];
        double pi = opi + alpha * m_Dir[entryParameter];

        double upr = m_pParamGroup->GetParamPtr(entryParameter)->GetUpperBoundTransformed();
        double lwr = m_pParamGroup->GetParamPtr(entryParameter)->GetLowerBoundTransformed();
        if (pi > upr) { alpha = (upr - opi) / m_Dir[entryParameter]; }
        if (pi < lwr) { alpha = (lwr - opi) / m_Dir[entryParameter]; }
    }

    return alpha;
}
