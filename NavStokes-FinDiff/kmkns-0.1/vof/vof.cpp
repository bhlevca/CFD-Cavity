#include <cmath>
#include <iostream>

#include "vof.h"

// If the slope of a line is less than this, it is rounded to zero.
const double SLOPE_EPSILON = 1e-8;
// If the difference between two values is less than this, they are "equal".
const double EQUAL_EPSILON = 1e-12;

inline void normalize(double& nx, double& ny)
{
	double l = sqrt(nx * nx + ny * ny);
	if(l != 0.0){nx /= l; ny /= l;}
}

const double PI = 3.1415926535897931;

namespace vof
{

// (-x,y)   (x,y)
//    +-\------+
//    |  \     |
//    |   \    |      x = dx/2
//    |    \   |      y = dy/2
//    |     \  |
//    +------\-+
// (-x,-y)   (x,-y)
// Calculate the fraction of the square area below the given line.
// The normal (nx, ny) is assumed to be normalized.
// 'd' is the signed distance from the origin to the line.
double fraction(double nx, double ny, double d, double dx, double dy)
{
	nx = vof::abs(nx);
	ny = vof::abs(ny);

	// Check if the line actually passes through the cell.
	// If top right corner is below the line, the entire cell is below the line.
	if(d + 0.5 * (nx * dx + ny * dy) <= 0.0)
		return 1.0;
	// If bottom left corner is above the line, the entire cell is above the line.
	if(d - 0.5 * (nx * dx + ny * dy) >= 0.0)
		return 0.0;

	// Treat vertical and horizontal lines as special cases.
	if(nx <= SLOPE_EPSILON)
		return (0.5 - d / dy); // needs verification
	if(ny <= SLOPE_EPSILON)
		return (0.5 - d / dx); // needs verification

	/*
	     Y
	     |\
	     | \
	     |  \
	     |   \
	     |    \
	     +-----T--+
	     |      \ |
	     |       \|
	     |        R
	     |        |\
	     +--------+-X
	*/

	// Calculate the line's intersection with the axes
	double x = (0.5 * dy * ny - d) / (dx * nx) + 0.5;
	double y = (0.5 * dx * nx - d) / (dy * ny) + 0.5;

	// Calculate the area of triangle (0, 0)-X-Y
	double totalArea = 0.5 * x * y;

	// Calculate relative area of top triangle (0, 1)-T-Y if Y is above (0, 1).
	double topFrac = 0.0;
	if(y > 1.0)
		topFrac = sqr((y-1.0)/y);

	// Calculate relative area of right triangle (1, 0)-X-R if R is to the right of (1, 0).
	double rightFrac = 0.0;
	if(x > 1.0)
		rightFrac = sqr((x-1.0)/x);

	return totalArea * (1.0 - topFrac - rightFrac);
}

/*
   (-x,y)   (x,y)
      +-\------+
      |  \     |
      |   \    |      x = dx/2
      |    \   |      y = dy/2
      |     \  |
      +------\-+
   (-x,-y)   (x,-y)
   Calculates the distance between the origin and the line whose normal is (nx, ny).
   The returned number will be negative if the origin is below the line, positive if
   above the line and zero if exactly on the line.
   The fraction of the rectangle area below the line will be equal to 'frac'.
   The normal (nx, ny) is assumed to be normalized.
*/
double distance(double nx, double ny, double frac, double dx, double dy)
{
	nx = vof::abs(nx);
	ny = vof::abs(ny);

	// Treat vertical and horizontal lines as special cases.
	if(nx <= SLOPE_EPSILON)
		return (0.5 - frac) * dy; // needs verification
	if(ny <= SLOPE_EPSILON)
		return (0.5 - frac) * dx; // needs verification

	// Calculate the signed distance from the centre of the square (0.5, 0.5)
	// to four potensial lines, each passing through one of the corners of the square.
	double dLL = 0.5 * (dx * nx + dy * ny); // LL = lower left corner
	double dUL = 0.5 * (dx * nx - dy * ny); // UL = upper left corner
	double dLR = -dUL; // LR = lower right corner
	double dUR = -dLL; // UR = upper right corner

	// Calculate the covered area if the line passes through the
	// upper left or lower right corner of the square.
	double fracUL = fraction(nx, ny, dUL, dx, dy);
	double fracLR = fraction(nx, ny, dLR, dx, dy);

	// If the area below both lines exceed the target fraction,
	// the final line must pass through the left and bottom edges of the square.
	if(frac <= fracUL && frac <= fracLR)
	{
		if(fracUL > fracLR)
		{
			swap(fracUL, fracLR);
			swap(dUL, dLR);
		}
		return dLL - sqrt(sqr(dUL - dLL) * frac / fracUL);
	}
	// If the area below both lines is less than the target fraction,
	// the final line must pass through the right and top edges of the square.
	else if(frac >= fracUL && frac >= fracLR)
	{
		if(fracLR < fracUL)
		{
			swap(fracUL, fracLR);
			swap(dUL, dLR);
		}
		return dUR + sqrt(sqr(dLR - dUR) * (1.0 - frac) / (1.0 - fracLR));
	}
	// Otherwise, the line must cross the square from left to right or top to bottom.
	else
	{
		// If the line passes through both corners, dUL==dLR, fracLR==fracUL==0.5.
		if(vof::abs(fracLR - fracUL) < EQUAL_EPSILON)
			return 0.5 * (dUL + dLR);

		// Linearly interpolate
		double theta = (frac - fracUL) / (fracLR - fracUL);
		return dLR * theta + dUL * (1.0 - theta);
	}
}

// Calculate the flux through the left edge of the cell.
// The flux equals the area below the line and within a distance of 'vdt' from
// the left edge divided by (dx * dy) to end up in the range [0, 1].
// Normal (nx, ny) is assumed to be normalized.
// 'd' is the signed distance from the centre of the square to the line.
// 'vdt' is "velocity times delta time" and must be positive and less than 'dx'.
double leftFlux(double nx, double ny, double d, double dx, double dy, double vdt)
{
	if(vdt == 0.0)
		return 0.0;
	return fraction(nx, ny, d + 0.5 * (vdt - dx) * nx, vdt, dy) * vdt / dx;
}

// etc.
double rightFlux(double nx, double ny, double d, double dx, double dy, double vdt)
{
	if(vdt == 0.0)
		return 0.0;
	return fraction(nx, ny, d - 0.5 * (vdt - dx) * nx, vdt, dy) * vdt / dx;
}

double topFlux(double nx, double ny, double d, double dx, double dy, double vdt)
{
	if(vdt == 0.0)
		return 0.0;
	return fraction(nx, ny, d - 0.5 * (vdt - dy) * ny, dx, vdt) * vdt / dy;
}

double bottomFlux(double nx, double ny, double d, double dx, double dy, double vdt)
{
	if(vdt == 0.0)
		return 0.0;
	return fraction(nx, ny, d + 0.5 * (vdt - dy) * ny, dx, vdt) * vdt / dy;
}

// Calculate the area of the intersection between the right half of the coordinate
// system (x > 0) and the circle centred at (d, 0) with radius 'r'.
double halfCircleIntersection(double d, double r)
{
	// solve equation y^2 + d^2 = r^2
	double r2d2 = r * r - d * d;
	if(r2d2 <= 0.0)
		return (d > 0.0 ? PI * r * r : 0.0);

	double y = sqrt(r * r - d * d); // positive and negative
	return r * r * (PI - atan2(y, d)) + y * d;
}

// Calculate the area of the intersection between the first quadrant (x, y > 0)
// and the circle centred at (cx, cy) with radius 'r'.
double quadrantCircleIntersection(double cx, double cy, double r)
{
	double area = 0.0;
	if(cx * cx + cy * cy >= r * r)
	{
		if(cy > 0.0) area += halfCircleIntersection(cx, r);
		if(cx > 0.0) area += halfCircleIntersection(cy, r);
		if(cx > 0.0 && cy > 0.0) area -= PI * r * r;
		return area;
	}
	double x = sqrt(r * r - cy * cy);
	double y = sqrt(r * r - cx * cx);
	double theta = acos(-(x * cx + cy * y) / (r * r));
	return 0.5 * (r * r * theta + (x + cx) * cy + (y + cy) * cx);
}

// Calculate the area of the intersection between the rectangle (x1, y1)-(x2, y1)-(x2, y2)-(x1, y2)
// and the circle centred at (cx, cy) with radius 'r'.
// 'x1' < 'x2', 'y1' < 'y2'
double rectCircleIntersection(double x1, double y1, double x2, double y2, double cx, double cy, double r)
{
	return quadrantCircleIntersection(cx - x1, cy - y1, r) +
		quadrantCircleIntersection(cx - x2, cy - y2, r) -
		quadrantCircleIntersection(cx - x1, cy - y2, r) -
		quadrantCircleIntersection(cx - x2, cy - y1, r);
}

// Calculate the length of the arc which is part of both the right half of the
// coordinate system (x > 0) and the circle centred at (d, 0) with radius 'r'.
double halfArcIntersection(double d, double r)
{
	// solve equation y^2 + d^2 = r^2
	double r2d2 = r * r - d * d;
	if(r2d2 <= 0.0)
		return (d > 0.0 ? 2.0 * PI * r : 0.0);

	double y = sqrt(r * r - d * d); // positive and negative
	return 2.0 * r * (PI - atan2(y, d));
}

// Calculate the length of the arc which is part of both the first quadrant (x, y > 0)
// and the circle centred at (cx, cy) with radius 'r'.
double quadrantArcIntersection(double cx, double cy, double r)
{
	double length = 0.0;
	if(cx * cx + cy * cy >= r * r)
	{
		if(cy > 0.0) length += halfArcIntersection(cx, r);
		if(cx > 0.0) length += halfArcIntersection(cy, r);
		if(cx > 0.0 && cy > 0.0) length -= 2.0 * PI * r;
		return length;
	}
	double x = sqrt(r * r - cy * cy);
	double y = sqrt(r * r - cx * cx);
	double theta = acos(-(x * cx + cy * y) / (r * r));
	return r * theta;
}

// Calculate the length of the arc which is part of both the rectangle
// (x1, y1)-(x2, y1)-(x2, y2)-(x1, y2) and the circle centred at (cx, cy) with radius 'r'.
// 'x1' < 'x2', 'y1' < 'y2'
double rectArcIntersection(double x1, double y1, double x2, double y2, double cx, double cy, double r)
{
	return quadrantArcIntersection(cx - x1, cy - y1, r) +
		quadrantArcIntersection(cx - x2, cy - y2, r) -
		quadrantArcIntersection(cx - x1, cy - y2, r) -
		quadrantArcIntersection(cx - x2, cy - y1, r);
}

// Calculate the area of the intersection between the rectangle (x11, y11)-(x12, y11)-(x12, y12)-(x11, y12)
// and the rectangle (x21, y21)-(x22, y21)-(x22, y22)-(x21, y22).
// 'x11' < 'x12', 'y11' < 'y12', 'x21' < 'x22', 'y21' < 'y22'
double rectRectIntersection(double x11, double y11, double x12, double y12,
							double x21, double y21, double x22, double y22)
{
	double dx = vof::max(vof::min(x12, x22) - vof::max(x11, x21), 0.0);
	double dy = vof::max(vof::min(y12, y22) - vof::max(y11, y21), 0.0);
	return dx * dy;
}

// Calculate the length of the line which is part of both the line (x1, y)-(x2, y)
// and the circle around the origin with radius 'r'.
// 'x1' <= 'x2'
double lineCircleIntersection(double x1, double x2, double y, double r)
{
	// solve x^2 + y^2 = r^2
	double xSqr = r * r - y * y;
	if(xSqr <= 0.0)
		return 0.0; // line does not intersect circle
	double x = sqrt(xSqr);
	return vof::max(0.0, vof::min(x, x2) - vof::max(-x, x1));
}

// Calculate the gradient of a scalar field using finite differences.
// The gradient is calculated for the cell corresponding to 'f[4]'.
// Scalar values are laid out like this:
//
// y      |--dx--|
// ^
// +------+------+------+
// |      |      |      |
// | f[6] | f[7] | f[8] |
// |      |      |      |
// +------+------+------+ -+-
// |      |      |      |  |
// | f[3] | f[4] | f[5] |  dy
// |      |      |      |  |
// +------+------+------+ -+-
// |      |      |      |
// | f[0] | f[1] | f[2] |
// |      |      |      |
// +------+------+------+--> x
void gradient(double& gx, double& gy, double* f, double dx, double dy)
{
	double r2p1 = 1.0 + sqr(dx / dy);
	double invr2p1 = 1.0 + sqr(dx / dy);
	double factor_x = 1.0/(dx * (4.0 + 2.0 * invr2p1));
	double factor_y = 1.0/(dy * (4.0 + 2.0 * r2p1));
	gx = factor_x * (f[8] + f[2] - f[6] - f[0] + invr2p1 * (f[5] - f[3]));
	gy = factor_y * (f[8] + f[6] - f[2] - f[0] + r2p1 * (f[7] - f[1]));
}

// Calculate the sum of squared errors caused by the given interface line across
// all nine cells. In the process, the line contant 'd' is also calculated and is returned.
double interfaceError(double& d, double nx, double ny, double* f, double dx, double dy)
{
	d = distance(nx, ny, f[4], dx, dy);
	double error = 0.0;
	error += sqr(f[0] - fraction(nx, ny, d - dx * nx - dy * ny, dx, dy));
	error += sqr(f[1] - fraction(nx, ny, d - dy * ny, dx, dy));
	error += sqr(f[2] - fraction(nx, ny, d + dx * nx - dy * ny, dx, dy));
	error += sqr(f[3] - fraction(nx, ny, d - dx * nx, dx, dy));
	error += sqr(f[5] - fraction(nx, ny, d + dx * nx, dx, dy));
	error += sqr(f[6] - fraction(nx, ny, d - dx * nx + dy * ny, dx, dy));
	error += sqr(f[7] - fraction(nx, ny, d + dy * ny, dx, dy));
	error += sqr(f[8] - fraction(nx, ny, d + dx * nx + dy * ny, dx, dy));
	return error;
}

// Returns true if the interface passes through the centre cell 'f[4]' or
// along one of its edges.
bool isInterfaceCell(double centre, double south, double west, double north, double east)
{
	return (centre > EQUAL_EPSILON && centre < 1.0 - EQUAL_EPSILON) ||
		(vof::abs(east - centre) >= 1.0 - 2.1 * EQUAL_EPSILON) ||
		(vof::abs(west - centre) >= 1.0 - 2.1 * EQUAL_EPSILON) ||
		(vof::abs(south - centre) >= 1.0 - 2.1 * EQUAL_EPSILON) ||
		(vof::abs(north - centre) >= 1.0 - 2.1 * EQUAL_EPSILON);
}

// Reconstruct the interface of volume-of-fluid field. This function will
// not necessarily return a normal parallel to the gradient returned by the gradient() function.
// The reconstruct() function will recreate a linear interface, while the gradient() function
// will not. The returned value is an L2-norm of the error due to the linear reconstruction, and
// will lie in the interval [0, 9].
// y      |--dx--|
// ^
// +------+------+------+
// |      |      |      |
// | f[6] | f[7] | f[8] |
// |      |      |      |
// +------+------+------+ -+-
// |      |      |      |  |
// | f[3] | f[4] | f[5] |  dy
// |      |      |      |  |
// +------+------+------+ -+-
// |      |      |      |
// | f[0] | f[1] | f[2] |
// |      |      |      |
// +------+------+------+--> x
double reconstruct(double& nx, double& ny, double& d, double* f, double dx, double dy)
{
	if(!isInterfaceCell(f[4], f[1], f[3], f[5], f[7]))
	{
		nx = ny = d = 0.0;
		return 0.0;
	}

	double c1 = f[0] + f[3] + f[6], c2 = f[1] + f[4] + f[7], c3 = f[2] + f[5] + f[8]; // column sums
	double r1 = f[0] + f[1] + f[2], r2 = f[3] + f[4] + f[5], r3 = f[6] + f[7] + f[8]; // row sums

	double bestError = 1e308;
	double nxs[6], nys[6], ds[6];

	// left difference
	nxs[0] = dy * (c1 - c2);
	nys[0] = (r1 > r3 ? dx : -dx);
	// right difference
	nxs[1] = dy * (c2 - c3);
	nys[1] = (r1 > r3 ? dx : -dx);
	// central difference
	nxs[2] = dy * (c1 - c3);
	nys[2] = 2.0 * (r1 > r3 ? dx : -dx);

	// bottom difference
	nys[3] = dx * (r1 - r2);
	nxs[3] = (c1 > c3 ? dy : -dy);
	// top difference
	nys[4] = dx * (r2 - r3);
	nxs[4] = (c1 > c3 ? dy : -dy);
	// central difference
	nys[5] = dx * (r1 - r3);
	nxs[5] = 2.0 * (c1 > c3 ? dy : -dy);

	for(int i = 0; i < 6; ++i)
	{
		normalize(nxs[i], nys[i]);
		// &f11 will point into the stack where f12, f13 etc. will follow.
		double error = interfaceError(ds[i], nxs[i], nys[i], f, dx, dy);
		if(error < bestError)
		{
			nx = nxs[i];
			ny = nys[i];
			d = ds[i];
			bestError = error;
		}
	}
	return bestError;
}

// Calculate the curvature (kappa=1/R) of an arc which passes through the three given points.
// The points should be given in clockwise order.
double curvatureFromPoints(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double ux = x2 - x1;
	double uy = y2 - y1;
	double vx = x3 - x2;
	double vy = y3 - y2;
	double wx = x1 - x3;
	double wy = y1 - y3;
	return 2.0 * (vx * uy - vy * ux) / sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy) * (wx * wx + wy * wy));
}

}
