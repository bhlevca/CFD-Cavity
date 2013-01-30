#ifndef VOF_H
#define VOF_H

namespace vof
{

const double CURVATURE_UNDEFINED = 1e308;

double fraction(double nx, double ny, double d, double dx, double dy);
double distance(double nx, double ny, double fraction, double dx, double dy);

double leftFlux(double nx, double ny, double d, double dx, double dy, double vdt);
double rightFlux(double nx, double ny, double d, double dx, double dy, double vdt);
double topFlux(double nx, double ny, double d, double dx, double dy, double vdt);
double bottomFlux(double nx, double ny, double d, double dx, double dy, double vdt);

double halfCircleIntersection(double d, double r);
double quadrantCircleIntersection(double cx, double cy, double r);
double rectCircleIntersection(double x1, double y1, double x2, double y2, double cx, double cy, double r);
double rectRectIntersection(double x11, double y11, double x12, double y12,
							double x21, double y21, double x22, double y22);
double halfArcIntersection(double d, double r);
double quadrantArcIntersection(double cx, double cy, double r);
double rectArcIntersection(double x1, double y1, double x2, double y2, double cx, double cy, double r);
double lineCircleIntersection(double x1, double x2, double y, double r);

void gradient(double& gx, double& gy, double* f, double dx, double dy);
double interfaceError(double& d, double nx, double ny, double* f, double dx, double dy);
double reconstruct(double& nx, double& ny, double& d, double* f, double dx, double dy);
double curvatureFromPoints(double x1, double y1, double x2, double y2, double x3, double y3);
bool isInterfaceCell(double centre, double south, double west, double north, double east);

inline double sqr(double x){return x * x;}
inline void swap(double& x, double& y){double t = x; x = y; y = t;}
inline double clamp(double x, double a, double b){return (x < a ? a : (x > b ? b : x));}
inline double max(double a, double b){return (a > b ? a : b);}
inline double min(double a, double b){return (a < b ? a : b);}
inline double sign(double x){return (x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0));}
inline double abs(double x){return (x > 0.0 ? x : -x);}

}

#endif
