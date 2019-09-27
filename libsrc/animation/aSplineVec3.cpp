#include "aSplineVec3.h"
#include <algorithm>
#include <Eigen/Dense>

#pragma warning(disable:4018)
#pragma warning(disable:4244)

ASplineVec3::ASplineVec3() : mInterpolator(new ALinearInterpolatorVec3())
{
}

ASplineVec3::~ASplineVec3()
{
  delete mInterpolator;
}

void ASplineVec3::setFramerate(double fps)
{
  mInterpolator->setFramerate(fps);
}

double ASplineVec3::getFramerate() const
{
  return mInterpolator->getFramerate();
}

void ASplineVec3::setLooping(bool loop)
{
  mLooping = loop;
}

bool ASplineVec3::getLooping() const
{
  return mLooping;
}

void ASplineVec3::setInterpolationType(ASplineVec3::InterpolationType type)
{
  double fps = getFramerate();

  delete mInterpolator;
  switch (type)
  {
  case LINEAR: mInterpolator = new ALinearInterpolatorVec3();
    break;
  case CUBIC_BERNSTEIN: mInterpolator = new ABernsteinInterpolatorVec3();
    break;
  case CUBIC_CASTELJAU: mInterpolator = new ACasteljauInterpolatorVec3();
    break;
  case CUBIC_MATRIX: mInterpolator = new AMatrixInterpolatorVec3();
    break;
  case CUBIC_HERMITE: mInterpolator = new AHermiteInterpolatorVec3();
    break;
  case CUBIC_BSPLINE: mInterpolator = new ABSplineInterpolatorVec3();
    break;
  };

  mInterpolator->setFramerate(fps);
  computeControlPoints();
  cacheCurve();
}

ASplineVec3::InterpolationType ASplineVec3::getInterpolationType() const
{
  return mInterpolator->getType();
}

void ASplineVec3::editKey(int keyID, const vec3& value)
{
  assert(keyID >= 0 && keyID < mKeys.size());
  mKeys[keyID].second = value;
  computeControlPoints();
  cacheCurve();
}

void ASplineVec3::editControlPoint(int ID, const vec3& value)
{
  assert(ID >= 0 && ID < mCtrlPoints.size()+2);
  if (ID == 0)
  {
    mStartPoint = value;
    computeControlPoints();
  }
  else if (ID == mCtrlPoints.size() + 1)
  {
    mEndPoint = value;
    computeControlPoints();
  }
  else mCtrlPoints[ID - 1] = value;
  cacheCurve();
}

void ASplineVec3::appendKey(double time, const vec3& value, bool updateCurve)
{
  mKeys.push_back(Key(time, value));

  if (mKeys.size() >= 2)
  {
    int totalPoints = mKeys.size();

    //If there are more than 1 interpolation point, set up the 2 end points to help determine the curve.
    //They lie on the tangent of the first and last interpolation points.
    vec3 tmp = mKeys[0].second - mKeys[1].second;
    double n = tmp.Length();
    mStartPoint = mKeys[0].second + (tmp / n) * n * 0.25;
    // distance to endpoint is 25% of distance between first 2 points

    tmp = mKeys[totalPoints - 1].second - mKeys[totalPoints - 2].second;
    n = tmp.Length();
    mEndPoint = mKeys[totalPoints - 1].second + (tmp / n) * n * 0.25;
  }

  if (updateCurve)
  {
    computeControlPoints();
    cacheCurve();
  }
}

void ASplineVec3::appendKey(const vec3& value, bool updateCurve)
{
  if (mKeys.size() == 0)
  {
    appendKey(0, value, updateCurve);
  }
  else
  {
    double lastT = mKeys[mKeys.size() - 1].first;
    appendKey(lastT + 1, value, updateCurve);
  }
}

void ASplineVec3::deleteKey(int keyID)
{
  assert(keyID >= 0 && keyID < mKeys.size());
  mKeys.erase(mKeys.begin() + keyID);
  computeControlPoints();
  cacheCurve();
}

vec3 ASplineVec3::getKey(int keyID)
{
  assert(keyID >= 0 && keyID < mKeys.size());
  return mKeys[keyID].second;
}

int ASplineVec3::getNumKeys() const
{
  return mKeys.size();
}

vec3 ASplineVec3::getControlPoint(int ID)
{
  assert(ID >= 0 && ID < mCtrlPoints.size()+2);
  if (ID == 0) return mStartPoint;
  else if (ID == mCtrlPoints.size() + 1) return mEndPoint;
  else return mCtrlPoints[ID - 1];
}

int ASplineVec3::getNumControlPoints() const
{
  return mCtrlPoints.size() + 2; // include endpoints
}

void ASplineVec3::clear()
{
  mKeys.clear();
}

double ASplineVec3::getDuration() const
{
  return mKeys[mKeys.size() - 1].first;
}

double ASplineVec3::getNormalizedTime(double t) const
{
  return (t / getDuration());
}

vec3 ASplineVec3::getValue(double t)
{
  if (mCachedCurve.size() == 0) return vec3();

  double dt = mInterpolator->getDeltaTime();
  int rawi = (int)(t / dt); // assumes uniform spacing
  int i = rawi % mCachedCurve.size();
  double frac = t - rawi * dt;
  int inext = i + 1;
  if (!mLooping) inext = std::min<int>(inext, mCachedCurve.size() - 1);
  else inext = inext % mCachedCurve.size();

  vec3 v1 = mCachedCurve[i];
  vec3 v2 = mCachedCurve[inext];
  vec3 v = v1 * (1 - frac) + v2 * frac;
  return v;
}

void ASplineVec3::cacheCurve()
{
  mInterpolator->interpolate(mKeys, mCtrlPoints, mCachedCurve);
}

void ASplineVec3::computeControlPoints()
{
  mInterpolator->computeControlPoints(mKeys, mCtrlPoints, mStartPoint, mEndPoint);
}

int ASplineVec3::getNumCurveSegments() const
{
  return mCachedCurve.size();
}

vec3 ASplineVec3::getCurvePoint(int i) const
{
  return mCachedCurve[i];
}

//---------------------------------------------------------------------
AInterpolatorVec3::AInterpolatorVec3(ASplineVec3::InterpolationType t) : mDt(1.0 / 120.0), mType(t)
{
}

void AInterpolatorVec3::setFramerate(double fps)
{
  mDt = 1.0 / fps;
}

double AInterpolatorVec3::getFramerate() const
{
  return 1.0 / mDt;
}

double AInterpolatorVec3::getDeltaTime() const
{
  return mDt;
}
void AInterpolatorVec3::interpolate(const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints, std::vector<vec3>& curve)
{
	vec3 val = 0.0;
	double u = 0.0;

	curve.clear();
	//system("cls");

	int numSegments = keys.size() - 1;
	for (int segment = 0; segment < numSegments; segment++)
	{ 
		for (double t = keys[segment].first; t < keys[segment + 1].first - FLT_EPSILON; t += mDt) // mDt is the timeStep ?
		{
			// TODO: Compute u, fraction of duration between segment and segmentnext, for example,
			// u = 0.0 when t = keys[segment-1].first  
			// u = 1.0 when t = keys[segment].first
			u = (t - keys[segment].first) / (keys[segment + 1].first - keys[segment].first);

			val = interpolateSegment(keys, ctrlPoints, segment, u);
			curve.push_back(val);
		}
	}
	// add last point
	if (keys.size() > 1)
	{
		u = 1.0;
		val = interpolateSegment(keys, ctrlPoints, numSegments - 1, u);
		curve.push_back(val);
	}
}

// u = [0.0, 1.0];
vec3 LerpHelper(vec3 iKey0Pos, vec3 iKey1Pos, double u)
{
	vec3 temp;
	for (int i = 0; i < 3; i++)
	{
		temp[i] = iKey0Pos[i] * (1 - u) + u * iKey1Pos[i];
	}
	return temp;
}

vec3 ALinearInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 curveValue(0, 0, 0);
	vec3 key0 = keys[segment].second;
	vec3 key1 = keys[segment + 1].second;

	// TODO: 
	//Step 1: Create a Lerp helper function
	//Step 2: Linear interpolate between key0 and key1 so that u = 0 returns key0 and u = 1 returns key1
	curveValue = LerpHelper(key0, key1, u);

	return curveValue;
}
vec3 ABernsteinInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 b0;
  vec3 b1;
  vec3 b2;
  vec3 b3;
  vec3 curveValue(0, 0, 0);
  // TODO: 
  // Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
  // Step2: Compute the interpolated value f(u) point using  Bernstein polynomials
  //system("CLS");
  //std::cout << "Index of the segment:" << segment << std::endl;
  // The segment index starts from 0;
  int j = segment;
  b0 = ctrlPoints[j * 4];
  b1 = ctrlPoints[j * 4 + 1];
  b2 = ctrlPoints[j * 4 + 2];
  b3 = ctrlPoints[j * 4 + 3];

  double B_0_3, B_1_3, B_2_3, B_3_3;

  B_0_3 = (1 - t) * (1 - t) * (1 - t);
  B_1_3 = 3 * t * (1 - t) * (1 - t);
  B_2_3 = 3 * t * t * (1 - t);
  B_3_3 = t * t * t;
  curveValue = b0 * B_0_3 + b1 * B_1_3 + b2 * B_2_3 + b3 * B_3_3;

  return curveValue;
}


vec3 ACasteljauInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 b0;
  vec3 b1;
  vec3 b2;
  vec3 b3;
  vec3 curveValue(0, 0, 0);

  // TODO: 
  // Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
  // Step2: Compute the interpolated value f(u) point using  deCsteljau alogithm
  int j = segment;
  b0 = ctrlPoints[j * 4];
  b1 = ctrlPoints[j * 4 + 1];
  b2 = ctrlPoints[j * 4 + 2];
  b3 = ctrlPoints[j * 4 + 3];

  vec3 b_0_1 = LerpHelper(b0, b1, 0.5);
  vec3 b_1_1 = LerpHelper(b1, b2, 0.5);
  vec3 b_2_1 = LerpHelper(b2, b3, 0.5);
  vec3 b_0_2 = LerpHelper(b_0_1, b_1_1, 0.5);
  vec3 b_1_2 = LerpHelper(b_1_1, b_2_1, 0.5);
  vec3 b_0_3 = LerpHelper(b_0_2, b_1_2, 0.5);

  double B_0_3, B_1_3, B_2_3, B_3_3;

  if (t < 0.5)
  {
	  double u = t / 0.5;
	  B_0_3 = (1 - u) * (1 - u) * (1 - u);
	  B_1_3 = 3 * u * (1 - u) * (1 - u);
	  B_2_3 = 3 * u * u * (1 - u);
	  B_3_3 = u * u * u;

	  curveValue = b0 * B_0_3 + b_0_1 * B_1_3 + b_0_2 * B_2_3 + b_0_3 * B_3_3;
  }
  else if (t >= 0.5) {
	  double v = (t - 0.5) / 0.5;
	  B_0_3 = (1 - v) * (1 - v) * (1 - v);
	  B_1_3 = 3 * v * (1 - v) * (1 - v);
	  B_2_3 = 3 * v * v * (1 - v);
	  B_3_3 = v * v * v;

	  curveValue = b_0_3 * B_0_3 + b_1_2 * B_1_3 + b_2_1 * B_2_3 + b3 * B_3_3;
  }

  return curveValue;
}

vec3 AMatrixInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 b0;
  vec3 b1;
  vec3 b2;
  vec3 b3;
  vec3 curveValue(0, 0, 0);

  // TODO: 
  // Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
  // Step2: Compute the interpolated value f(u) point using  matrix method f(u) = GMU
  // Hint: Using Eigen::MatrixXd data representations for a matrix operations
  int j = segment;
  b0 = ctrlPoints[j * 4];
  b1 = ctrlPoints[j * 4 + 1];
  b2 = ctrlPoints[j * 4 + 2];
  b3 = ctrlPoints[j * 4 + 3];

  Eigen::Matrix<Eigen::Vector3d, -1, 1> gVector(4);
  Eigen::MatrixXd mMatrix(4, 4);
  //Eigen::Matrix<Eigen::Vector3d, 1, -1> gVector(4);
  Eigen::Vector4d uVector(1, t, t * t, t * t * t);

  gVector(0) = Eigen::Vector3d(b0[0], b0[1], b0[2]);
  gVector(1) = Eigen::Vector3d(b1[0], b1[1], b1[2]);
  gVector(2) = Eigen::Vector3d(b2[0], b2[1], b2[2]);
  gVector(3) = Eigen::Vector3d(b3[0], b3[1], b3[2]);

  mMatrix << 1, -3, 3, -1,
			 0, 3, -6, 3,
			 0, 0, 3, -3,
			 0, 0, 0, 1;
  
  Eigen::Vector4d temp = mMatrix * uVector;
  for (int i = 0; i < 4; i++)
  {
	  Eigen::Vector3d ele = temp(i) * gVector(i);
	  curveValue += vec3(ele(0), ele(1), ele(2));
  }

  return curveValue;
}

vec3 AHermiteInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 p0 = keys[segment].second;
  vec3 p1 = keys[segment + 1].second;
  vec3 q0 = ctrlPoints[segment]; // slope at p0
  vec3 q1 = ctrlPoints[segment + 1]; // slope at p1
  vec3 curveValue(0, 0, 0);

  // TODO: Compute the interpolated value h(u) using a cubic Hermite polynomial  
  double t_3 = t * t * t;
  double t_2 = t * t;

  double b1 = 2 * t_3 - 3 * t_2 + 1;
  double b2 = -2 * t_3 + 3 * t_2;
  double b3 = t_3 - 2 * t_2 + t;
  double b4 = t_3 - t_2;

  curveValue = b1 * p0 + b2 * p1 + b3 * q0 + b4 * q1;
  return curveValue;
}




// Used to calculate the value of Nj_n, which also is a recurisive function.
double Nj_n(std::vector<double> knotVec, int j, int n, double t)
{
	// Find almbda j and almbda j+1.
	double almbda_j = knotVec[j];
	double almbda_jadd1 = knotVec[j + 1];
	double almbda_jaddnadd1 = knotVec[j + n + 1];
	double almbda_jaddn = knotVec[j + n];

	// Return 0.0, if the ti is not located at the range of [almbda_j, almbda_j+n+1].
	// Check whether the input t is a ti.
	bool isTi = false;
	for (int i = 0; i < knotVec.size(); i++)
	{
		double dis = t - knotVec[i];
		if (abs(dis) < DBL_EPSILON)
		{
			isTi = true;
			break;
		}
	}

	if (isTi)
	{
		if ((t < almbda_j || t > almbda_jaddnadd1))
		{
			return 0.0;
		}
	}

	// Return 0.0 or 1.0 when n = 0.
	if (n == 0)
	{
		if (t < almbda_jadd1 && t >= almbda_j)
		{
			return 1.0;
		}
		else
		{
			return 0.0;
		}
	}

	double currentNj_n = (t - almbda_j) * Nj_n(knotVec, j, n - 1, t) / (almbda_jaddn - almbda_j) + (almbda_jaddnadd1 - t) * Nj_n(knotVec, j + 1, n - 1, t) / (almbda_jaddnadd1 - almbda_jadd1);

	return currentNj_n;
}

double Nj_n2(std::vector<double> knotVec, int j, int n, double t)
{
	// Check whether the input t is a ti.
	bool isTi = false;
	for (int i = 0; i < knotVec.size(); i++)
	{
		double dis = t - knotVec[i];
		if (abs(dis) < DBL_EPSILON)
		{
			isTi = true;
			break;
		}
	}

	if ((j + n + 1) >= knotVec.size())
	{
		return 0.0;
	}

	// Find almbda j and almbda j+1.
	double almbda_j = knotVec[j];
	double almbda_jadd1 = knotVec[j + 1];
	double almbda_jaddnadd1 = knotVec[j + n + 1];
	double almbda_jaddn = knotVec[j + n];

	// Return 0.0, if the ti is not located at the range of [almbda_j, almbda_j+n+1].
	

	if (isTi)
	{
		if ((t < almbda_j || t > almbda_jaddnadd1))
		{
			return 0.0;
		}
	}

	// Return 0.0 or 1.0 when n = 0.
	if (n == 0)
	{
		if (t < almbda_jadd1 && t >= almbda_j)
		{
			return 1.0;
		}
		else
		{
			return 0.0;
		}
	}

	double currentNj_n = (t - almbda_j) * Nj_n(knotVec, j, n - 1, t) / (almbda_jaddn - almbda_j) + (almbda_jaddnadd1 - t) * Nj_n(knotVec, j + 1, n - 1, t) / (almbda_jaddnadd1 - almbda_jadd1);

	return currentNj_n;
}

// Used to calculate the derivative of N.
double dN(std::vector<double> knotVec, int j, int n, double t, int l)
{
	if (l == 0)
	{
		double check = Nj_n(knotVec, j, n, t);
		return check;
	}
	double almbda_j = knotVec[j];
	double almbda_jadd1 = knotVec[j + 1];
	double almbda_jaddnadd1 = knotVec[j + n + 1];
	double almbda_jaddn = knotVec[j + n];

	// Return 0.0 or 1.0 when n = 0.
	if (n == 0)
	{
		if (t < almbda_jadd1 && t >= almbda_j)
		{
			return 1.0;
		}
		else
		{
			return 0.0;
		}
	}

	double Ele1 = dN(knotVec, j, n - 1, t, l - 1) / (almbda_jaddn - almbda_j);
	double Ele2 = dN(knotVec, j + 1, n - 1, t, l - 1) / (almbda_jaddnadd1 - almbda_jadd1);

	double current_dN = n * (Ele1 - Ele2);
	return current_dN;
}


vec3 ABSplineInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double t)
{
  vec3 curveValue(0, 0, 0);

  // Hint: Create a recursive helper function N(knots,n,j,t) to calculate BSpline basis function values at t, where
  //     knots = knot array
  //	   n = degree of the spline curves (n =3 for cubic)
  //     j = curve interval on knot vector in which to interpolate
  //     t = time value	

  // Step 1: determine the index j
  // Construct Knot Vector;
  std::vector<double> knotVec;
  double tInterval = keys[1].first - keys[0].first; // It should be 1.0;
  int degree = 3;
  for (int i = 0; i < degree; i++)
  {
	  int coeff = i - degree;
	  knotVec.push_back(keys[0].first + coeff * tInterval);
  }

  double knotVecEle = 0.0;
  for (int i = 0; i < keys.size(); i++, knotVecEle += 1.0)
  {
	  knotVec.push_back(knotVecEle);
  }

  for (int i = 0; i < degree; i++)
  {
	  int coeff = 1 + i;
	  knotVec.push_back(keys[keys.size() - 1].first + coeff * tInterval);
  }

  int j = 0;
  for (;; j++)
  {
	  if (t >= knotVec[j] && t < knotVec[j + 1])
	  {
		  break;
	  }
  }

  // Step 2: compute the n nonzero Bspline Basis functions N given j
  // Step 3: get the corresponding control points from the ctrlPoints vector
  // Step 4: compute the Bspline curveValue at time t
  
  for (int k = 0; k < 4; k++)
  {
	  //std::cout << "t:" << t << std::endl;
	  //std::cout << "j:" << j << std::endl;
	  //std::cout << "j - k:" << j - k << std::endl;
	  //double currentN = Nj_n2(knotVec, j - k, 3, 1);
	  //double currentN = Nj_n2(knotVec, 4, 3, 1.0);
	  double currentN = Nj_n2(knotVec, j - k, 3, t);
	  //std::cout << "t:" << t << std::endl;
	  //std::cout << "j:" << j << std::endl;
	  //std::cout << "j - k:" << j - k << std::endl;
	  /*vec3 currentCtrPoint = ctrlPoints[j - k];
	  curveValue[0] += currentCtrPoint[0] * currentN;
	  curveValue[1] += currentCtrPoint[1] * currentN;
	  curveValue[2] += currentCtrPoint[2] * currentN;*/
  }
  
  return curveValue;
}

void ACubicInterpolatorVec3::computeControlPoints(
  const std::vector<ASplineVec3::Key>& keys,
  std::vector<vec3>& ctrlPoints,
  vec3& startPoint, vec3& endPoint)
{
  ctrlPoints.clear();
  if (keys.size() <= 1) return;

  // 这个 i = 1 有点玄学 ==> 似乎是为了让这个i在含义上能够表示segment标号？ #0 segment, #1 segment, .... #N-1 segment；
  // 这个循环执行输入点数量少一次也是有点意思 ==> 相当于每一个segment执行一次循环；
  for (int i = 1; i < keys.size(); i++)
  {
    vec3 b0, b1, b2, b3;

    // TODO: compute b0, b1, b2, b3
	
	int j = i - 1;
	vec3 s0, s1;

	b0 = keys[j].second;
	b3 = keys[j + 1].second;

	if ((j - 1) < 0) {
		s0 = keys[1].second - keys[0].second;
	}
	else {
		s0 = (keys[j + 1].second - keys[j - 1].second) / 2;
	}

	if ((j + 2) > (keys.size() - 1)) {
		s1 = keys[j + 1].second - keys[j].second;
	}
	else {
		s1 = (keys[j + 2].second - keys[j].second) / 2;
	}

	b2 = b3 - s1 / 3;
	b1 = b0 + s0 / 3;
	
	ctrlPoints.push_back(b0);
	ctrlPoints.push_back(b1);
	ctrlPoints.push_back(b2);
	ctrlPoints.push_back(b3);
  }
}

void AHermiteInterpolatorVec3::computeControlPoints(
  const std::vector<ASplineVec3::Key>& keys,
  std::vector<vec3>& ctrlPoints,
  vec3& startPoint, vec3& endPoint)
{
  ctrlPoints.clear();
  if (keys.size() <= 1) return;

  int numKeys = keys.size();

  // TODO: 
  // For each key point pi, compute the corresonding value of the slope pi_prime.
  // Hints: Using Eigen::MatrixXd for a matrix data structures, 
  // this can be accomplished by solving the system of equations AC=D for C.
  // Don't forget to save the values computed for C in ctrlPoints
  // For clamped endpoint conditions, set 1st derivative at first and last points (p0 and pm) to s0 and s1, respectively
  // For natural endpoints, set 2nd derivative at first and last points (p0 and pm) equal to 0

  // Step 1: Initialize A
  // Step 2: Initialize D
  // Step 3: Solve AC=D for C
  // Step 4: Save control points in ctrlPoints

  // Note: be careful to the conversion between vec3 and Eigen::Vector3xxxxx.
  // Node: the default values in Eigen::Matrix / Vector are extreme large negative values rather than 0. Thus, you have to init them.
  
  Eigen::Matrix<Eigen::Vector3d, -1, 1> dVector(numKeys);
  Eigen::MatrixXd aMatrix(numKeys, numKeys);


  // Push every elements into the dArray.
  vec3 temp = 3 * (keys[1].second - keys[0].second);
  dVector(0)(0) = temp[0];
  dVector(0)(1) = temp[1];
  dVector(0)(2) = temp[2];
  
  temp = 3 * (keys[numKeys - 1].second - keys[numKeys - 2].second);
  dVector(numKeys - 1)(0) = temp[0];
  dVector(numKeys - 1)(1) = temp[1];
  dVector(numKeys - 1)(2) = temp[2];

  for (int i = 2; i < numKeys; i++)
  {
	  temp = 3 * (keys[i].second - keys[i - 2].second);
	  dVector(i - 1)(0) = temp[0];
	  dVector(i - 1)(1) = temp[1];
	  dVector(i - 1)(2) = temp[2];
  }

  // Push every elements into the aMatrix.
  // Set the first row.
  for (int i = 0; i < numKeys; i++)
  {
	  if (i == 0)
	  {
		  aMatrix(0, 0) = 2;
	  }
	  else if (i == 1)
	  {
		  aMatrix(0, 1) = 1;
	  }
	  else {
		  aMatrix(0, i) = 0;
	  }
  }
  // Set the last row.
  for (int i = 0; i < numKeys; i++)
  {
	  if (i == (numKeys - 2))
	  {
		  aMatrix(numKeys - 1, numKeys - 2) = 1;
	  }
	  else if (i == (numKeys - 1)) {
		  aMatrix(numKeys - 1, numKeys - 1) = 2;
	  }
	  else {
		  aMatrix(numKeys - 1, i) = 0;
	  }
  }
  // Set the rows at the middle.
  for (int i = 1; i < numKeys - 1; i++)
  {
	  for (int j = 0; j < numKeys; j++)
	  {
		  if (j != (i - 1) && j != (i) && j != (i + 1))
		  {
			  aMatrix(i, j) = 0;
		  }
		  else {
			  aMatrix(i, i - 1) = 1;
			  aMatrix(i, i) = 4;
			  aMatrix(i, i + 1) = 1;
		  }
	  }
  }
  
  // Compute the C = (A^(-1))D.
  // Init cVector.
  Eigen::Matrix<Eigen::Vector3d, -1, 1> cVector(numKeys);
  Eigen::Vector3d initEigenVec3(0, 0, 0);
  for (int i = 0; i < numKeys; i++)
  {
	  cVector(i) = initEigenVec3;
  }

  Eigen::MatrixXd aMatrixInverse(numKeys, numKeys);
  aMatrixInverse = aMatrix.inverse();

  for (int i = 0; i < numKeys; i++)
  {
	  Eigen::Vector3d tempCElement;
	  
	  // Extract #i row from the inverse aMatrix.
	  Eigen::VectorXd tempRow(numKeys);
	  for (int j = 0; j < numKeys; j++)
	  {
		  tempRow(j) = aMatrixInverse(i, j);
	  }

	  // Compute the #i row of the cVector.
	  // Use #j col of the #i row of the inverse matrix multiplies the #j row of the dVector.
	  // A.k.a use #j element of the tempRow multiplies the #j row of the dVector.
	  Eigen::Matrix<Eigen::Vector3d, -1, 1> middleValuesVector(numKeys);
	  for (int j = 0; j < numKeys; j++)
	  {
		  middleValuesVector(j) = tempRow(j) * dVector(j);
	  }
	  for (int j = 0; j < numKeys; j++)
	  {
		  vec3 check = vec3(middleValuesVector(j)(0), middleValuesVector(j)(1), middleValuesVector(j)(2));
		  cVector(i) = cVector(i) + middleValuesVector(j);
	  }
  }
  
  // Push all elements in the cVector into the ctrlPoints vector.
  
  for (int i = 0; i < numKeys; i++)
  {
	  vec3 check = vec3(vec3(cVector(i)(0), cVector(i)(1), cVector(i)(2)));
	  ctrlPoints.push_back(check);
  }
}

// Used to calculate the basis spline N
void ABSplineInterpolatorVec3::computeControlPoints(
  const std::vector<ASplineVec3::Key>& keys,
  std::vector<vec3>& ctrlPoints,
  vec3& startPt, vec3& endPt)
{
  ctrlPoints.clear();
  if (keys.size() <= 1) return;

  // TODO: c
  // Hints: 
  // 1. use Eigen::MatrixXd to calculate the control points by solving the system of equations AC=D for C

  // 2. Create a recursive helper function dN(knots,n,t,l) to calculate derivative BSpline values at t, where
  //     knots = knot array
  //	   n = degree of the spline curves (n =3 for cubic)
  //     j = interval on knot vector in which to interpolate
  //     t = time value
  //     l = derivative (l = 1 => 1st derivative)

  // Step 1: Calculate knot vector using a uniform BSpline
  //         (assume knots are evenly spaced 1 apart and the start knot is at time = 0.0)

  // 1.1: Construct the knotVec, and lambdaVec.

  std::vector<double> knotVec;
  double tInterval = keys[1].first - keys[0].first; // It should be 1.0;
  int degree = 3;
  for (int i = 0; i < degree; i++)
  {
	  int coeff = i - degree;
	  knotVec.push_back(keys[0].first + coeff * tInterval);
  }

  double knotVecEle = 0.0;
  for (int i = 0; i < keys.size(); i++, knotVecEle += 1.0)
  {
	  knotVec.push_back(knotVecEle);
  }

  for (int i = 0; i < degree; i++)
  {
	  int coeff = 1 + i;
	  knotVec.push_back(keys[keys.size() - 1].first + coeff * tInterval);
  }
  /**/
  //std::cout << knotVec[0] << std::endl;

  // Step 2: Calculate A matrix  for a natural BSpline
  //         (Set 2nd derivative at t0 and tm to zero, where tm is the last point knot; m = #segments)
  
  Eigen::MatrixXd aMatrix(keys.size() + 2, keys.size() + 2);
  // Init the Matrix.
  for (int i = 0; i < keys.size() + 2; i++)
  {
	  for (int j = 0; j < keys.size() + 2; j++)
	  {
		  aMatrix(i, j) = 0.0;
	  }
  }
  
  // Set the value
  // i is the row index, j is the column index.
  for (int i = 0; i < keys.size() + 2; i++)
  {
	  if (i == 0)
	  {
		  // The first row.
		  for (int j = 0; j < 4; j++)
		  {
			  aMatrix(i, j) = dN(knotVec, j, 3, keys[0].first, 2);
		  }
	  }
	  else if(i == keys.size() + 1)
	  {
		  // The last row.
		  for (int j = keys.size() - 2; j < keys.size() + 2; j++)
		  {
			  aMatrix(i, j) = dN(knotVec, j, 3, keys[keys.size() - 1].first, 2);
		  }
	  }
	  else {
		  // The rows at the middle.
		  for (int j = 0; j < keys.size() + 2; j++)
		  {
			  aMatrix(i, j) = Nj_n(knotVec, j, 3, keys[i - 1].first);
		  }
	  }
  }/**/
  //double check = Nj_n(knotVec, 0, 3, keys[0].first);
  //std::cout << check << std::endl;
  // double dN(std::vector<double> knotVec,int j, int n,double t,int l)
  //double check = dN(knotVec, 2, 3, keys[0].first, 2);
  //std::cout << check << std::endl;

  // Step 3: Calculate  D matrix composed of our target points to interpolate
  // Init the D vector.
  
  Eigen::Matrix<Eigen::Vector3d, -1, 1> dVector(keys.size() + 2);
  Eigen::Vector3d initEigenVec3(0, 0, 0);
  for (int i = 0; i < keys.size() + 2; i++)
  {
	  dVector(i) = initEigenVec3;
  }

  // Put key points into it.
  for (int i = 1; i < keys.size() + 1; i++)
  {
	  dVector(i) = Eigen::Vector3d(keys[i - 1].second[0], keys[i - 1].second[1], keys[i - 1].second[2]);
  }
  /**/
  // Step 4: Solve AC=D for C 
  // Get the A inverse.
  
  Eigen::MatrixXd aMatrixInverse(keys.size() + 2, keys.size() + 2);
  aMatrixInverse = aMatrix.inverse();
  /**/
  // Compute the C = A^(-1)D.
  
  Eigen::Matrix<Eigen::Vector3d, -1, 1> cVector(keys.size() + 2);
  for (int i = 0; i < keys.size() + 2; i++)
  {
	  cVector(i) = initEigenVec3;
  }

  for (int i = 0; i < keys.size() + 2; i++)
  {
	  Eigen::Vector3d tempDElement;

	  // Extract #i row from the inverse aMatrix.
	  Eigen::VectorXd tempRow(keys.size() + 2);
	  for (int j = 0; j < keys.size() + 2; j++)
	  {
		  tempRow(j) = aMatrixInverse(i, j);
	  }

	  // Compute the #i row of the cVector.
	  // Use #j col of the #i row of the inverse matrix multiplies the #j row of the dVector.
	  // A.k.a use #j element of the tempRow multiplies the #j row of the dVector.
	  Eigen::Matrix<Eigen::Vector3d, -1, 1> middleValuesVector(keys.size() + 2);
	  for (int j = 0; j < keys.size() + 2; j++)
	  {
		  middleValuesVector(j) = tempRow(j) * dVector(j);
	  }
	  for (int j = 0; j < keys.size() + 2; j++)
	  {
		  vec3 check = vec3(middleValuesVector(j)(0), middleValuesVector(j)(1), middleValuesVector(j)(2));
		  cVector(i) = cVector(i) + middleValuesVector(j);
	  }
  }
  /**/

  // Step 5: save control points in ctrlPoints
  for (int i = 0; i < keys.size() + 2; i++)
  {
	  vec3 check = vec3(vec3(cVector(i)(0), cVector(i)(1), cVector(i)(2)));
	  ctrlPoints.push_back(check);
  }/**/
}
