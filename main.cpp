#include <stdc++.h>
#define eps 1e-7
#define inf 1e300
#define PI 3.14159625358979323l
using namespace std;

// TODO: eliminate the dependency of the algorithm on non-vertical lines. 
// Also, personally I think this needed to be refactored later:
// The code is not concise enough in my opinion.
// But for now at least it works.

// Some functions of triangles to be added: altitudes (and etc.), excenter, (mixtillinear circle), similar triangles, etc.

// Graphics. How to generate graphics on points & lines & circles? At least support
// -> Scaling...
// How about SFML / others? (I just think SFML can be used immediately)

// BTW. How to name the points? or triangles?
// Maybe this can considered as a combination of name tag and fundamental structures of points, circles, lines, etc.

// Parsing of 'problems'

struct point {
	long double x, y;

	point() {};
	point(long double _x, long double _y) : x(_x), y(_y) {};

	point operator+		(const point& o) const { return { x + o.x, y + o.y }; }
	point operator+=	(const point& o) { x += o.x, y += o.y; }
	point operator-		(const point& o) const { return { x - o.x, y - o.y }; }
	point operator-=	(const point& o) { x -= o.x, y -= o.y; }
	point operator*		(long double n) const { return { x * n, y * n }; }
	friend point operator* (long double n, const point& A) { return A * n; }
	point operator*=	(long double n) { x *= n, y *= n; }
	point operator/		(long double n) const { return { x / n, y / n }; }
	point operator/=	(long double n) { x /= n, y /= n; }
	bool operator== (point o) { return (x - o.x) < eps && (y - o.y) < eps; }

	friend point sectional_formula(point A, point B, long double s, long double t) { // s : t, both s and t can be negative
		return ((A * t) + (B * s)) / (s + t);
	}

	friend point mid_point(point A, point B) {
		return (A + B) / 2;
	}

	void print() const {
		cout << '(' << x << ", " << y << ')' << endl;
	}

	friend std::ostream &operator<<(std::ostream &os, point const &P) {
		return os << '(' << P.x << ", " << P.y << ')';
	}
};

struct Ratio {
	long double p = 0, q = 1;
	Ratio(long double _p, long double _q) : p(_p), q(_q) { long double g = max(p, q); p /= g; q /= g; if (p < 0) { p = -p, q = -q; } };
	Ratio() {};
	bool operator== (const Ratio& r) const {
		return (p * r.q - q * r.p) < eps;
	}

	void print() const {
		cout << "Ratio: " << p << ":" << q << endl;
	}

	friend std::ostream &operator<<(std::ostream &os, Ratio const &ratio) {
		return os << "Ratio: " << ratio.p << ":" << ratio.q << endl;
	}
};

long double dis(point A, point B) {
	return sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}

long double dis_square(point A, point B) {
	return (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y);
}

struct line {
	long double a, b, c; // ax + by + c = 0
	long double l = -inf, r = inf; // for line segments & rays. Crucial in checking intersections.
	line(point A, point B) : a(A.y - B.y), b(B.x - A.x), c(A.x * B.y - A.y * B.x) {};
	line(point A, long double m) : a(m), b(-1), c(A.y - m * A.x) {};

	long double slope() { return -a / b; }
	bool on(point p) {
		return abs(a * p.x + b * p.y + c) <= eps;
	}

	long double sub(point p) {
		return a * p.x + b * p.y + c;
	}

	point unit(long double l) { // vector of length l, with the direction of the line
		return { l * cos(atan(slope())), l * sin(atan(slope())) };
	}

	void print() const {
		cout << "Equation of line: " << a << "x + " << b << "y + " << c << " = 0" << endl;
	}

	friend std::ostream &operator<<(std::ostream &os, line const &l) {
		return os << "Equation of line: " << l.a << "x + " << l.b << "y + " << l.c << " = 0" << endl;
	}
};

long double dis(line l, point A) {
	return abs(l.sub(A)) / sqrt(l.a * l.a + l.b * l.b);
}

bool on(line l, point A) {
	return dis(l, A) <= eps;
}

struct lineSegment : public line {
	lineSegment(point A, point B) : line(A, B) {
		l = A.x, r = B.x;
	}
};

struct ray : public line {
	long double angle;
	ray(point A, point B) : line(A, B) { // A->B-> 
		l = A.x, r = inf;
		angle = atan2(B.y - A.y, B.x - A.x);
	}
};

point intersection(line l_1, line l_2, int cnt = 0) {
	long double x = (l_1.b * l_2.c - l_2.b * l_1.c) / (l_2.b * l_1.a - l_1.b * l_2.a);
	long double y = (l_1.a * l_2.c - l_2.a * l_1.c) / (l_2.a * l_1.b - l_1.a * l_2.b); 
	return {x, y};
}

line perpendicular_bisector(point A, point B) {
	point M = (A + B) / 2;
	long double m = (long double) -1 / line(A, B).slope();
	return line(M, m);
}

line perpendicular(point A, line l) {
	return line(A, -1 / l.slope());
}

line internal_angle_bisector(point A, point B, point C) {
	point C_prime = sectional_formula(B, C, dis(A, B), dis(B, C) - dis(A, B));
	return perpendicular_bisector(A, C_prime);
}

line external_angle_bisector(point A, point B, point C) {
	return perpendicular(B, internal_angle_bisector(A, B, C));
}

struct triangle {
	point A, B, C;
	triangle() {};
	triangle(point A, point B, point C) : A(A), B(B), C(C) {};

	triangle(line a, line b, line c) {
		A = intersection(b, c);
		B = intersection(c, a);
		C = intersection(a, b);
	}

	long double s() {
		return (dis(A, B) + dis(B, C) + dis(C, A)) / 2;
	}

	point orthocentre() {
		return intersection(perpendicular(A, line(B, C)), perpendicular(B, line(C, A)));
	}

	point incentre() {
		return intersection(internal_angle_bisector(A, B, C), internal_angle_bisector(B, C, A));
	}

	point centroid() {
		return (A + B + C) / 3;
	}

	point circumcentre() {
		return intersection(perpendicular_bisector(A, B), perpendicular_bisector(B, C));
	}

	long double area() { // Heron's formula
		long double a = dis(B, C), b = dis(C, A), c = dis(A, B);
		long double s = (a + b + c) / 2;
		return sqrt(s * (s - a) * (s - b) * (s - c));
	}
};

struct circle {
	point centre;
	long double radius;
	circle(point O, long double rad) { centre = O; radius = rad; }

	circle(point A, point B, point C) {
		point O = triangle(A, B, C).circumcentre();
		long double rad = dis(O, A);
		circle(O, rad);
	}

	circle(triangle Triangle) {
		point O = Triangle.circumcentre();
		long double rad = dis(O, Triangle.A);
		centre = O, radius = rad;
	}

	line tangent(point P) {
		return line(P, -1 / line(centre, P).slope());
	}

	long double power(point X) {
		return (dis(centre, X) * dis(centre, X)) - radius * radius;
	}

	bool on(point X) {
		return abs(power(X)) < eps;
	}
};

struct directed_angle {
	long double rad = 0;
	void scale() { while (rad > PI) rad -= PI; while (rad < -PI) rad += PI; }
	directed_angle() {};
	directed_angle(point A, point B, point C) { rad = ray(B, C).angle - ray(B, A).angle; scale(); };
};

vector<point> intersection(circle w, line l, int cnt = 0) {
	vector<point> solutions;
	if (w.radius <= dis(l, w.centre) - eps) {
		; // 0 solution
	} else if (abs(w.radius - dis(l, w.centre)) < eps) {
		solutions = { intersection(l, perpendicular(w.centre, l)) }; // 1 solution
	} else {
		point M = intersection(l, perpendicular(w.centre, l)); // 2 solution
		long double r = sqrt(w.radius * w.radius - dis(l, w.centre) * dis(l, w.centre));
		solutions.push_back(M + l.unit(r));
		solutions.push_back(M - l.unit(r));
	}
	return solutions;
}

circle circumcircle(triangle Triangle) {
	return circle(Triangle);
}

circle incircle(triangle Triangle) {
	point I = Triangle.incentre();
	long double r = Triangle.area() / Triangle.s();
	return circle(I, r);
}

point mid_point_of_arc(circle w, point A, point B, bool cnt = 0) { // cnt = 0 for minor arc; cnt = 1 for major arc
	vector<point> vp = intersection(w, internal_angle_bisector(A, w.centre, B));
	if (vp.size() <= 1) {
		cerr << "Error: Points are not on the circle" << endl;
		return { inf, inf };
	}else if(vp.size() == 2){
		if (dis(A, vp[0]) < dis(B, vp[1]) ^ cnt) {
			return vp[0];
		} else {
			return vp[1];
		}
	}
}

bool collinear(point A, point B, point C) {
	return line(A, B).on(C);
}

bool parallel(line l_1, line l_2) {
	return abs(l_1.slope() - l_2.slope()) <= eps;
}

bool concyclic(point A, point B, point C, point D) {
	if (collinear(A, B, C)) return false;
	return circle(A, B, C).on(D);
}

line radical_axis(circle w_1, circle w_2) {
	point P = sectional_formula(w_1.centre, w_2.centre, (w_1.radius * w_1.radius - w_2.radius * w_2.radius) / 2 / dis(w_1.centre, w_2.centre) + dis(w_1.centre, w_2.centre) / 2, (w_2.radius * w_2.radius - w_1.radius * w_1.radius) / 2 / dis(w_1.centre, w_2.centre) + dis(w_1.centre, w_2.centre) / 2);
	long double m = -1 / line(w_1.centre, w_2.centre).slope();
	return line(P, m);
}

point radical_centre(circle w_1, circle w_2, circle w_3) {
	return intersection(radical_axis(w_1, w_2), radical_axis(w_2, w_3));
}

vector<point> intersection(circle w_1, circle w_2, int cnt = 0) {
	return intersection(w_1, radical_axis(w_1, w_2), cnt);
}

struct vector2d {
	long double x, y;

	vector2d() {};
	vector2d(long double _x, long double _y) : x(_x), y(_y) {};

	vector2d operator+		(const vector2d& o) const { return { x + o.x, y + o.y }; }
	vector2d operator+=		(const vector2d& o) { x += o.x, y += o.y; }
	vector2d operator-		(const vector2d& o) const { return { x - o.x, y - o.y }; }
	vector2d operator-=		(const vector2d& o) { x -= o.x, y -= o.y; }
	vector2d operator*		(long double n) const { return { x * n, y * n }; }
	friend vector2d operator* (long double n, const vector2d& A) { return A * n; }
	vector2d operator*=		(long double n) { x *= n, y *= n; }
	vector2d operator/		(long double n) const { return { x / n, y / n }; }
	vector2d operator/=		(long double n) { x /= n, y /= n; }
	friend long double cross(vector2d A, vector2d B) { return A.x * B.y - A.y * B.x; }

};

int32_t main() {

	point A(5, 3), B(1, -1), C(3, 6);
	triangle T(A, B, C);
	cout << "H" << T.orthocentre() << endl;
	cout << "G" << T.centroid() << endl;
	cout << "I" << T.incentre() << endl;
	cout << "O" << T.circumcentre() << endl;
	cout << "HG:GO = " << Ratio(dis(T.orthocentre(), T.centroid()), dis(T.centroid(), T.circumcentre())) << endl;
	cout << A + B << endl;
	cout << A - B << endl;
	cout << A / 2 << endl;
	cout << B / 2 << endl;
	cout << dis(A, B) << endl;
	cout << dis_square(A, B) << endl;
	cout << sectional_formula(A, B, 1, 3) << endl;
	cout << sectional_formula(A, B, -1, 3) << endl;
}
