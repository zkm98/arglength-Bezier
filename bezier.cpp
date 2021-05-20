#include<iostream>
#include<time.h>
#include<vector>
#include<queue>
#include<unordered_map>
#include<string>

using namespace std;

struct NumericalIntegration
{
	template <typename F>
	static double
		simpson_3_8(F&& derivative_func, const double& L, const double& R)
	{
		double mid_L = (2 * L + R) / 3.0, mid_R = (L + 2 * R) / 3.0;
		return (derivative_func(L) +
			3.0 * derivative_func(mid_L) +
			3.0 * derivative_func(mid_R) +
			derivative_func(R)) * (R - L) / 8.0;
	}

	template <typename F>
	static double
		adaptive_simpson_3_8(F&& derivative_func,
			const double& L, const double& R, const double& eps = 0.0001)
	{
		const double mid = (L + R) / 2.0;
		double ST = simpson_3_8(derivative_func, L, R),
			SL = simpson_3_8(derivative_func, L, mid),
			SR = simpson_3_8(derivative_func, mid, R);
		double ans = SL + SR - ST;
		if (fabs(ans) <= 15.0 * eps)  return SL + SR + ans / 15.0;
		return adaptive_simpson_3_8(derivative_func, L, mid, eps / 2.0) +
			adaptive_simpson_3_8(derivative_func, mid, R, eps / 2.0);
	}
};
struct Point {
	double x;
	double y;
	double norm() {
		return sqrt(x * x + y * y);
	}
	Point(double x_, double y_){
		x = x_;
		y = y_;
	}
	Point(){
		x = 0.0;
		y = 0.0;
	}
	Point operator*=(float a) {
		x *= a;
		y *= a;
		return *this;
	}
	Point operator* (double a) {
		double xx = x * a;
		double yy = y * a;
		return Point(xx, yy);
	}
	Point operator+(const Point a) {
		double xs = a.x + x;
		double ys = a.y + y;
		return Point(xs, ys);
	}
	Point operator-(const Point a) {
		double xs = -a.x + x;
		double ys = -a.y + y;
		return Point(xs, ys);
	}
	Point operator=(const Point& a) {
		x = a.x;
		y = a.y;
		return *this;
	}
	Point operator+=(const Point& a) {
		x += a.x;
		y += a.y;
		return *this;
	}
};

class BezierCurve2arglenth {
public:
	Point p0, p1, p2;
	BezierCurve2arglenth(Point a, Point b, Point c) {
		p0 = a;
		p1 = b;
		p2 = c;
	}
	BezierCurve2arglenth() {}
	Point at(double t) {
		return p0 * (pow(1 - t, 2)) + p1 * 2 * t * (1 - t) + p2 * t * t;
	}
	Point dev(double t) {
		double k1 = 2 * t;
		double k2 = 2 - k1;
		return ((p2 - p1) * k1 + (p1 - p0) * k2);
	}
	Point dev2(double t) {
		return (p2 + p0 - p1 * 2.0) * 2.0;
	}
	double total_length() {
		auto df = [&](double t) -> double
		{
			return this->dev(t).norm();
		};
		return NumericalIntegration::adaptive_simpson_3_8(df, 0, 1);
	}
	inline double simple_common(double x,double c,double d) {
		double k1 = abs((x + c) / d);
		double s = sqrt(k1 * k1 + 1);
		return k1 * s + log(k1 + s);
	}
	double length_with_t(double t) {
		double k1 = p1.x - p0.x;
		double k2 = p1.y - p1.y;
		double b1 = p0.x + p2.x - 2 * p1.x;
		double b2 = p0.y + p2.y - 2 * p2.y;

		double c = (k1 * b1 + k2 * b2) / (b1 * b1 + b2 * b2);
		double d = (k1 * b2 - k2 * b1) / (b1 * b1 + b2 * b2);
		double q1 = sqrt(b1 * b1 + b2 * b2) * abs(d) / 2;
		//cout << "when t equals " << t << " the common is " << simple_common(t,c,d) << endl;
		return (simple_common(t, c, d) - simple_common(0, c, d))*q1;
	}
	vector<Point> getPoints(int k) {
		vector<Point> ans(k + 2);
		ans[0] = p0;
		ans[k + 1] = p2;
		if (k < 1)return ans;

		vector<double> dist(k + 2, 0.0);
		vector<double> t_array(k + 2, 0.0);
		for (int i = 1; i < k + 2; i++) {
			dist[i] = (double(i)) / (double)(k + 1);
			t_array[i] = (double(i)) / (double)(k + 1);
			ans[i] = at(dist[i]);
		}
		double avg_distance = total_length() / ((double)(k + 1));

		int max_iter_time = 10000;
		for (int i = 0; i < max_iter_time; i++) {
			for (int j = 1; j <= k; ++j)dist[j] = (ans[j] - ans[j - 1]).norm();
			double offset = 0;
			for (int j = 1; j <= k; ++j) {
				const double error_dis = dist[j] - avg_distance;
				offset += error_dis;

				double first_dev = dev(t_array[j]).norm();
				double second_dev = dev2(t_array[j]).norm();
				double numerator = offset * first_dev;
				double denominator = offset * second_dev + first_dev * first_dev;

				t_array[j] -= numerator / denominator;
				ans[j] = at(t_array[j]);
			}
		}
		return ans;
	}
};

class BezierCurve {
public:
	vector<int> jiecheng;
	vector<Point> pts;
	int n;
	BezierCurve(vector<Point> pts_) {
		pts = pts_;
		n = pts.size();
		jiecheng = vector<int>(n + 1, 1);
		for (int i = 2; i <= n; i++) {
			jiecheng[i] = i * jiecheng[i - 1];
		}
	}
	Point at(double t) {
		Point ans;
		for (int i = 0; i < pts.size(); i++) {
			ans += pts[i] * mypow(1 - t, n-1 - i) * mypow(t, i) * comb(n-1, i);
		}
		return ans;
	}
	double total_length() {
		auto df = [&](double t) -> double
		{
			return this->dev(t).norm();
		};
		return NumericalIntegration::adaptive_simpson_3_8(df, 0, 1);
	}
	int comb(int n, int r) {
		if (r > n || r < 0)return 0;
		return jiecheng[n] / jiecheng[r] / jiecheng[n - r];
	}
	double mypow(double r, int k,int dev = 0) {
		if (k > n || k < -dev)return 0;
		return pow(r, k);
	}
	Point dev(double t) {
		int N = n - 1;
		if (t == 0)return pts[0] * N;
		if (t == 1.0)return pts[n - 1] * N;
		Point ans;
		for (int i = 0; i < n; i++) {
			ans = ans + pts[i] * comb(N, i) * (i * mypow(1 - t, N - i) * mypow(t, i - 1) - (double)(N - i) * mypow(1 - t, N - i - 1) * mypow(t, i));
		}
		return ans;
	}
	Point dev2(double t) {
		Point ans;
		int N = n - 1;

		for (int i = 0; i < n; ++i) {
			ans = ans + pts[i] * comb(N - 1, i) * (	(N - i) * (N - i - 1) * mypow(1 - t, N - i - 2)	- 2 * (N - i) * i * mypow(1 - t, N - i - 1) * mypow(t, i - 1)+ i * (N - 1) * mypow(1 - t, N - i) * mypow(t, i - 2));
		}
		return ans;
	}

	vector<Point> get_arg_points(int k, double iter_eps = 1e-5) {
		vector<Point> ans(k+2);
		ans[0] = pts[0];
		ans[k + 1] = pts.back();
		if (k < 1) {
			return ans;
		}
		vector<double> dist(k + 2, 0.0);
		vector<double> t_array(k + 2, 0.0);
		for (int i = 1; i < k + 2; i++) {
			dist[i] = (double(i)) / (double)(k+1);
			t_array[i] = (double(i)) / (double)(k + 1);
			ans[i] = at(dist[i]);
		}
		double avg_distance = total_length() / ((double)(k + 1));

		int max_iter_time = 10000;
		for (int i = 0; i < max_iter_time; i++) {
			for (int j = 1; j <= k; ++j)dist[j] = (ans[j] - ans[j - 1]).norm();
			double offset = 0;
			double abs_err = 0.0;
			for (int j = 1; j <= k; ++j) {
				const double error_dis = dist[j] - avg_distance;
				offset += error_dis;
				abs_err = abs(error_dis);
				double first_dev = dev(t_array[j]).norm();
				double second_dev = dev2(t_array[j]).norm();
				double numerator = offset * first_dev;
				double denominator = offset * second_dev + first_dev * first_dev;

				t_array[j] -= numerator / denominator;
				ans[j] = at(t_array[j]);
			}
			if (abs_err / k <= iter_eps)break;
		}
		return ans;
	}

	vector<Point> get_points(int k) {
		double delta_t = 1.0 / (k + 1);
		vector<Point> ans(k + 2);
		ans[0] = pts[0];
		for (int i = 1; i < k + 2; i++)ans[i] = at(delta_t * i);
		return ans;
	}
};

class Bezier_smoothL3 {
public:

	Bezier_smoothL3() {
	}

	Point at(vector<Point>ct, double t) {
		Point ans = ct[0] * mypow(1 - t, 3) + (ct[1] * (1 - t) + ct[2] * t) * (1 - t) * t + ct[3] * mypow(t, 3);
		return ans;
	}
	double mypow(double r, int k) {
		if (k < 0)return 0;
		return pow(r, k);
	}
	Point dev(vector<Point>&ct, double t) {
		Point ans = (ct[1] - ct[0]) * (3 * mypow(1 - t, 2)) + (ct[2] - ct[1]) * (6 * t * (1 - t)) + (ct[3] - ct[2]) * (3 * t * t);
		return ans;
	}
	Point dev2(vector<Point>& ct, double t) {
		Point ans;
		ans = (ct[0] + ct[2] - ct[1] * 2) * (6 * (1 - t)) + (ct[1] + ct[3] - ct[2] * 2) * (6 * t);
		return ans;
	}
	vector<Point> fit_curve(vector<Point>& pts, const double& maxE) {
		int sz = pts.size();
		vector<Point> ans;

		Point left_tan = pts[1] - pts[0];
		left_tan = left_tan * (1 / left_tan.norm());
		Point right_tan = pts[sz - 2] - pts[sz - 1];
		right_tan = right_tan * (1 / right_tan.norm());
		return fit_cubic(pts, 0, sz-1, left_tan, right_tan, maxE);
	}

	vector<Point> fit_cubic(vector<Point>& pts, int first, int last, Point left_tangent, Point right_tangent, const double& maxErr) {
		vector<Point> ans;
		int npts = last - first + 1;
		if (npts == 2) {
			double dis = (pts[last] - pts[first]).norm() / 3;
			ans.push_back(pts[first]);
			ans.push_back(pts[first] + left_tangent * dis);
			ans.push_back(pts[last] + right_tangent * dis);
			ans.push_back(pts[last]);
			return ans;
		}
		vector<double> u = chordLengthParameterize(pts,first,last);
		vector<Point> bezeCurv = generateBezier(pts, first,last, u, left_tangent, right_tangent);
		double maxError;
		int splitPoint;
		computeMaxError(pts, first, last , bezeCurv, u, maxError, splitPoint);
		printf("max error: %f", maxError);
		if (maxError < maxErr)return bezeCurv;
		if (maxErr > 1 && maxError < maxErr * maxErr) {
			for (int j = 0; j < 20; ++j) {
				vector<double> up = reparameterize(bezeCurv, pts, u, first, last);
				bezeCurv = generateBezier(pts, first, last, up, left_tangent, right_tangent);
				computeMaxError(pts, first, last,bezeCurv, u, maxError, splitPoint);
				if (maxError < maxErr)
					return bezeCurv;
				u = up;
			}
		}
		vector<Point> beziers;
 		Point t = pts[splitPoint - 1] - pts[splitPoint + 1];
		Point centerTangent = t * (1 / t.norm());
		vector<Point> ptsF(pts.begin(), pts.begin() + splitPoint);
		vector<Point> ptsB(pts.begin() + splitPoint, pts.end());

		beziers = fit_cubic(ptsF,first, splitPoint,left_tangent,centerTangent,maxErr);
		vector<Point> k = fit_cubic(ptsB, splitPoint, last,centerTangent, right_tangent, maxErr);
		beziers.insert(beziers.end(), k.begin(), k.end());
		return beziers;
	}

	vector<Point> generateBezier(vector<Point>& pts, const int&first, const int last, vector<double> parameters, Point left_tag, Point right_tag) {
		vector<Point> ans(4);
		ans[0] = pts[first]; ans[3] = pts[last];

		//if (pts.size() != parameters.size())cerr << "This is wrong between the parameter size and points size in function generateBezier!\n";
		int sz = last - first + 1;


		vector<vector<vector<double>>> A(sz, vector<vector<double>>(2, vector<double>(2, 0.0)));
		for (int i = 0; i < sz; ++i) {
			double u = parameters[i];
			double u_ = 1 - u;
			Point temp = left_tag * (3 * u * u_ * u_);
			A[i][0][0] = temp.x;
			A[i][0][1] = temp.y;
			temp = right_tag * (3 * u * u * u_);
			A[i][1][0] = temp.x;
			A[i][1][1] = temp.y;
		}

		vector<vector<double>> C(2, vector<double>(2, 0));
		vector<double> X(2, 0.0);

		for (int i = 0; i < sz; ++i) {
			C[0][0] += A[i][0][0] * A[i][0][0] + A[i][0][1] * A[i][0][1];
			C[0][1] += A[i][0][0] * A[i][1][0] + A[i][0][1] * A[i][1][1];
			C[1][0] += A[i][0][0] * A[i][1][0] + A[i][0][1] * A[i][1][1];
			C[1][1] += A[i][1][0] * A[i][1][0] + A[i][1][1] * A[i][1][1];

			Point temp = pts[i+first] - at({ pts[first],pts[first],pts[last],pts[last] }, parameters[i]);

			X[0] += A[i][0][0] * temp.x + A[i][0][1] * temp.y;
			X[1] += A[i][1][0] * temp.x + A[i][1][1] * temp.y;
		}

		double det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
		double det_C0_X = C[0][0] * X[1] - C[1][0] * X[0];
		double det_X_C1 = X[0] * C[1][1] - X[1] * C[0][1];
		
		double alpha_l = det_C0_C1 == 0 ? 0.0 : det_X_C1 / det_C0_C1;
		double alpha_r = det_C0_C1 == 0 ? 0.0 : det_C0_X / det_C0_C1;

		double segLength = (pts[first] - pts[last]).norm();
		double eps = segLength * 1e-6;
		if (alpha_l < eps || alpha_r < eps) {
			ans[1] = ans[0] + left_tag * (segLength / 3.0);
			ans[2] = ans[3] + right_tag * (segLength / 3.0);
		}
		else {
			ans[1] = ans[0] + left_tag * alpha_l;
			ans[2] = ans[3] + right_tag * alpha_r;
		}
		return ans;
	}

	vector<double> chordLengthParameterize(vector<Point>& pts, const int& first, const int &last) {
		vector<double> us = { 0 };
		for (int i = first+1; i <= last; ++i) {
			us.push_back(us.back() + (pts[i] - pts[i - 1]).norm());
		}
		for (int i = 0; i < us.size(); ++i) {
			us[i] /= us.back();
		}
		return us;
	}

	vector<double> reparameterize(vector<Point>& ctrls, vector<Point>& pts, vector<double> parameters,const int&first, const int&last) {
		vector<double> ans;
		//if (pts.size() != parameters.size())cerr << "The size of points and parameters is not Equall in reparameter!";
		int sz = last - first + 1;
		for (int i = 0; i < sz; ++i) {
			ans.push_back(newtonRaphsonRootFind(ctrls, pts[i+first], parameters[i]));
		}
		return ans;
	}

	double newtonRaphsonRootFind(vector<Point>& ctrls,const Point& p, const double& u) {
		Point d = at(ctrls, u) - p;
		Point dev1 = dev(ctrls, u);
		Point dev22 = dev2(ctrls, u);
		double numerator = d.x * dev1.x + d.y * dev1.y;
		double denominator = dev1.x * dev1.x + dev1.y * dev1.y + d.x * dev22.x + d.y * dev22.y;
		if (numerator == 0.0)return u;
		return u - numerator / denominator;
	}

	void computeMaxError(vector<Point>& pts, const int& first, const int & last, vector<Point>& ctrls, vector<double> params, double& maxDis, int& splitPoint) {
		double md = 0.0;
		int sp = (first + last) / 2;
		int sz = last - first +1;
		for (int i = 0; i < sz; ++i) {
			double dist = (at(ctrls, params[i])-pts[i+first]).norm();
			if (dist * dist > md) {
				md = dist * dist;
				sp = i;
			}
		}
		maxDis = md;
		splitPoint = sp;
	}

};


void main(){
	Bezier_smoothL3 b;
	vector<Point> pts = {
		{ 0.0, 0.0 },
		{ 0.0, 0.5 },
		{ 1.1, 1.4 },
		{ 2.1, 1.6 },
		{ 3.2, 1.1 },
		{ 4.0, 0.2 },
		{ 4.0, 0.0 }
	};
	cout << "pts: \n";
	for (int i = 0; i < pts.size(); ++i) {
		cout << "[" << pts[i].x << "," << pts[i].y << "],";
	}
	cout << "\n";
	double error = 4.0;
	vector<Point> ans = b.fit_curve(pts, error);
	for (int i = 0; i < ans.size() / 4; ++i) {
		vector<Point> pts;
		for (int j = 0; j < 4; ++j) {
			cout << "[" << ans[i * 4 + j].x << "," << ans[i * 4 + j].y << "],";
			pts.push_back(ans[i * 4 + j]);
		}
		cout << "\n";
		BezierCurve bc(pts);
		vector<Point> as = bc.get_points(40);
		for (int j = 0; j < 40; ++j) {
			cout << "[" << as[j].x << "," << as[j].y << "],";
		}
		cout << endl;
	}
}