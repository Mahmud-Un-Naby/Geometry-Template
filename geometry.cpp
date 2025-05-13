#include <bits/stdc++.h>
#define enl '\n'
#define all(x) begin(x), end(x)
using namespace std;
using i64 = long long;

struct Point;
struct Line;
struct LineSegment;
struct Polygon;

struct Point {
    long long x, y;
    Point();
    Point(long long xx, long long yy);
    friend istream& operator>>(istream& cin, Point& p);
    friend ostream& operator<<(ostream& cout, Point const& p);
    Point operator-(Point const& other) const;
    void operator-=(Point const& other);
    long long operator*(Point const& other) const;
    long long dot(Point const& other) const;
    long long triangle2(Point const& b, Point const& c) const;
    long double triangle(Point const& b, Point const& c);
    int orientation(Point const& b, Point const& c) const;
    long double angle(Point const& b, Point const& c) const;
    long double distanceToLine(LineSegment const&) const;
    long double distanceToLine(Line const&) const;
    long double distanceToLineSegment(LineSegment const&) const;
    long double distanceToRay(LineSegment const&) const;
};

struct Line {
    long long A, B, C; // coefficient of Ax + By + C = 0
    friend istream& operator>>(istream& cin, Line& L);
    bool isParallel(Line const& other) const;
    Point intersection(Line const& other) const;
};

struct LineSegment {
    Point A, B;
    LineSegment();
    friend istream& operator>>(istream& cin, LineSegment& ls);
    long double distanceToLineSegment(LineSegment const&) const;
    bool withinRectangle(Point const& p) const;
    bool intersects(LineSegment const& other) const;
};

struct Polygon {
    vector<Point> a;
    long double Area;
    Polygon(int n);
    friend istream& operator>>(istream& cin, Polygon& p);
    friend ostream& operator<<(ostream& cout, Polygon& p);
    bool on(Point const& p);
    bool in(Point const& p);
    long double area();
};

// Member method definitions

// Point methods
Point::Point() : x(0), y(0) { }
Point::Point(long long xx, long long yy) : x(xx), y(yy) { }
istream& operator>>(istream& cin, Point& p) {
    cin >> p.x >> p.y;
    return cin;
}
ostream& operator<<(ostream& cout, Point const& p) {
    cout << p.x << " " << p.y;
    return cout;
}
Point Point::operator-(Point const& other) const {
    return Point{ x - other.x, y - other.y };
}
void Point::operator-=(Point const& other) {
    x -= other.x;
    y -= other.y;
}
long long Point::operator*(Point const& other) const {
    return x * other.y - y * other.x;
}
long long Point::dot(Point const& other) const {
    return x * other.x + y * other.y;
}
long long Point::triangle2(Point const& b, Point const& c) const {
    return (b - *this) * (c - *this);
}
long double Point::triangle(Point const& b, Point const& c) {
    return abs(this->triangle2(b, c) * 0.5);
}
int Point::orientation(Point const& b, Point const& c) const {
    auto res = this->triangle2(b, c);
    if (res > 0) return 2;
    if (res < 0) return 1;
    return 0;
}
long double Point::angle(Point const& b, Point const& c) const {
    Point ab = b - *this;
    Point ac = c - *this;
    long long dot_product = ab.dot(ac);
    long double magnitude_ab = sqrtl(ab.x * ab.x + ab.y * ab.y);
    long double magnitude_ac = sqrtl(ac.x * ac.x + ac.y * ac.y);
    if (magnitude_ab == 0 || magnitude_ac == 0) {
        return 0;
    }
    return acos(dot_product / (magnitude_ab * magnitude_ac));
}
long double Point::distanceToLine(LineSegment const& ls) const {
    Point AB = ls.B - ls.A;
    if (AB.x == 0 && AB.y == 0) {
        // Segment is just a point
        return sqrtl(((*this) - ls.A).dot((*this) - ls.A));
    }
    Point AP = *this - ls.A;
    long long cross = AB * AP;
    long double base = sqrtl(AB.dot(AB));
    return abs(cross) / base;
}
long double Point::distanceToLine(Line const& L) const {
    return abs(L.A * x + L.B * y + L.C) / sqrtl(L.A * L.A + L.B * L.B);
}
long double Point::distanceToRay(LineSegment const& ray) const {
    Point AB = ray.B - ray.A;
    Point AP = *this - ray.A;
    long long dot1 = AP.dot(AB);
    if (dot1 < 0) {
        return sqrtl(AP.dot(AP));
    }
    long long cross = AB * AP;
    long double base = sqrtl(AB.dot(AB));
    return abs(cross) / base;
}
long double Point::distanceToLineSegment(LineSegment const& ls) const {
    Point AB = ls.B - ls.A;
    if (AB.x == 0 && AB.y == 0) {
        // Segment is a point
        return sqrtl(((*this) - ls.A).dot((*this) - ls.A));
    }
    Point AP = *this - ls.A;
    long long dot1 = AP.dot(AB);
    long long dot2 = AB.dot(AB);
    if (dot1 > 0 && dot1 < dot2) {
        // Projection falls strictly between A and B
        long long cross = AB * AP;
        long double base = sqrtl(dot2);
        return abs(cross) / base;
    }
    // Otherwise, return the distance to the nearest endpoint
    long double d1 = sqrtl(((*this) - ls.A).dot((*this) - ls.A));
    long double d2 = sqrtl(((*this) - ls.B).dot((*this) - ls.B));
    return min(d1, d2);
}

// Line methods
istream& operator>>(istream& cin, Line& L) {
    cin >> L.A >> L.B >> L.C;
    return cin;
}
bool Line::isParallel(Line const& other) const {
    return A * other.B == B * other.A;
}
Point Line::intersection(Line const& other) const {
    long long det = A * other.B - B * other.A;
    if (det == 0) return Point(-1e18, -1e18); // Lines are parallel or coincident
    long long x = (B * other.C - C * other.B) / det;
    long long y = (C * other.A - A * other.C) / det;
    return Point(x, y);
}

// LineSegment methods
LineSegment::LineSegment() : A(), B() { }
istream& operator>>(istream& cin, LineSegment& ls) {
    cin >> ls.A >> ls.B;
    return cin;
}
long double LineSegment::distanceToLineSegment(LineSegment const& ls) const {
    if (this->intersects(ls))
        return 0;
    long double d1 = ls.A.distanceToLineSegment(*this);
    long double d2 = ls.B.distanceToLineSegment(*this);
    long double d3 = A.distanceToLineSegment(ls);
    long double d4 = B.distanceToLineSegment(ls);
    return min({ d1, d2, d3, d4 });
}
bool LineSegment::withinRectangle(Point const& p) const {
    return min(A.x, B.x) <= p.x && p.x <= max(A.x, B.x) &&
        min(A.y, B.y) <= p.y && p.y <= max(A.y, B.y);
}
bool LineSegment::intersects(LineSegment const& other) const {
    int o1 = A.orientation(B, other.A);
    int o2 = A.orientation(B, other.B);
    int o3 = other.A.orientation(other.B, A);
    int o4 = other.A.orientation(other.B, B);
    if (o1 != o2 && o3 != o4)
        return true;
    if (o1 == 0 && withinRectangle(other.A)) return true;
    if (o2 == 0 && withinRectangle(other.B)) return true;
    if (o3 == 0 && other.withinRectangle(A)) return true;
    if (o4 == 0 && other.withinRectangle(B)) return true;
    return false;
}

// Polygon methods
Polygon::Polygon(int n) {
    a.resize(n);
    Area = -1;
}
istream& operator>>(istream& cin, Polygon& p) {
    for (auto& i : p.a) cin >> i;
    return cin;
}
ostream& operator<<(ostream& cout, Polygon& p) {
    for (auto& i : p.a) cout << i << enl;
    return cout;
}
bool Polygon::on(Point const& p) {
    int n = a.size();
    for (int i = 0; i < n; i++) {
        Point A = a[i], B = a[(i + 1) % n];
        if (p.x < min(A.x, B.x) || p.x > max(A.x, B.x))
            continue;
        if (p.y < min(A.y, B.y) || p.y > max(A.y, B.y))
            continue;
        if (p.orientation(A, B) == 0)
            return true;
    }
    return false;
}
bool Polygon::in(Point const& p) {
    if (on(p)) return false;
    int n = a.size();
    int cnt = 0;
    for (int i = 0; i < n; i++) {
        Point A = a[i], B = a[(i + 1) % n];
        if (A.y > B.y) swap(A, B);
        if (p.y > A.y && p.y <= B.y && (B - A) * (p - A) > 0) {
            cnt++;
        }
    }
    return cnt & 1;
}
long double Polygon::area() {
    if (Area >= 0)
        return Area;
    Area = 0;
    int n = a.size();
    for (int i = 2; i < n; i++) {
        Area += a[0].triangle2(a[i - 1], a[i]);
    }
    Area = abs(Area) * 0.5;
    return Area;
}

// Main function
void work() {
    Line a, b;
    cin >> a >> b;
    cout << a.intersection(b) << enl;
}

int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    while (t--) {
        work();
    }
}