// Copyright 2017 Nicolas Mellado
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// -------------------------------------------------------------------------- //
//
// Authors: Dror Aiger, Yoni Weill, Nicolas Mellado
//
// This file is part of the implementation of the 4-points Congruent Sets (4PCS)
// algorithm presented in:
//
// 4-points Congruent Sets for Robust Surface Registration
// Dror Aiger, Niloy J. Mitra, Daniel Cohen-Or
// ACM SIGGRAPH 2008 and ACM Transaction of Graphics.
//
// Given two sets of points in 3-space, P and Q, the algorithm applies RANSAC
// in roughly O(n^2) time instead of O(n^3) for standard RANSAC, using an
// efficient method based on invariants, to find the set of all 4-points in Q
// that can be matched by rigid transformation to a given set of 4-points in P
// called a base. This avoids the need to examine all sets of 3-points in Q
// against any base of 3-points in P as in standard RANSAC.
// The algorithm can use colors and normals to speed-up the matching
// and to improve the quality. It can be easily extended to affine/similarity
// transformation but then the speed-up is smaller because of the large number
// of congruent sets. The algorithm can also limit the range of transformations
// when the application knows something on the initial pose but this is not
// necessary in general (though can speed the runtime significantly).

// Home page of the 4PCS project (containing the paper, presentations and a
// demo): http://graphics.stanford.edu/~niloy/research/fpcs/fpcs_sig_08.html
// Use google search on "4-points congruent sets" to see many related papers
// and applications.
#define SUPER4PCS_USE_OPENMP
//#pragma comment(lib,"mclmcrrt.lib")
//#pragma comment(lib,"libmx.lib")
//#pragma comment(lib,"libmat.lib")
#ifndef _SUPER4PCS_ALGO_MATCH_4PCS_BASE_IMPL_
#define _SUPER4PCS_ALGO_MATCH_4PCS_BASE_IMPL_

#ifndef _SUPER4PCS_ALGO_MATCH_4PCS_BASE_
#include "super4pcs/algorithms/match4pcsBase.h"
#endif
#define OPENCV_TRAITS_ENABLE_DEPRECATED
#include <chrono>
#include <atomic>
#include <opencv2/opencv.hpp>
#include <time.h>

#include <vector>
#include <iostream>
#include <E:\Program Files\PCL 1.8.1\3rdParty\Eigen\eigen3\Eigen/Geometry>                 // MatrixBase.homogeneous()
#include <E:\Program Files\PCL 1.8.1\3rdParty\Eigen\eigen3\Eigen/SVD>                      // Transform.computeRotationScaling()
using namespace cv;
using namespace std;


namespace GlobalRegistration{

// The main 4PCS function. Computes the best rigid transformation and transfoms
// Q toward P by this transformation.
template <typename Sampler, typename Visitor>
Match4PCSBase::Scalar

Match4PCSBase::ComputeTransformation(const std::vector<Point3D>& P,
	                                       std::vector<Point3D::VectorType>& normals1,
	                                       std::vector<Point3D>* Q,
	                                      std::vector<Point3D::VectorType>* normals2,
                                     Eigen::Ref<MatrixType> transformation,
                                     const Sampler& sampler,
                                     const Visitor& v) {

  if (Q == nullptr) return kLargeNumber;
  if (P.empty() || Q->empty()) return kLargeNumber;

  init(P, normals1, *Q, *normals2, sampler);

  if (best_LCP_ != Scalar(1.))
    Perform_N_steps(number_of_trials_, transformation, Q, v);

#ifdef TEST_GLOBAL_TIMINGS
  Log<LogLevel::Verbose>( "----------- Timings (msec) -------------" );
  Log<LogLevel::Verbose>( " Total computation time  : ", totalTime   );
  Log<LogLevel::Verbose>( " Total verify time       : ", verifyTime  );
  Log<LogLevel::Verbose>( "    Kdtree query         : ", kdTreeTime  );
  Log<LogLevel::Verbose>( "----------------------------------------" );
#endif

  return best_LCP_;
}

//void linear_search(Mat& data, Mat& point, Mat& indices, Mat& dists, const int k) {
//	cv::flann::Index flannIndex(data, flann::LinearIndexParams(), cvflann::FLANN_DIST_L2);
//	cout << "Begin LinearSearch: " << endl;
//	clock_t tt1 = clock();
//	flannIndex.knnSearch(point, indices, dists, k, flann::SearchParams(64));
//	clock_t tt2 = clock();
//	cout << "Finish linerSearch: " << (tt2 - tt1)/ CLOCKS_PER_SEC << "s" << endl;
//}

template <typename Sampler>
void Match4PCSBase::init(
	const std::vector<Point3D>& P,
	const std::vector<Point3D::VectorType>& normals1,
	const std::vector<Point3D>& Q,
	const std::vector<Point3D::VectorType>& normals2,
	const Sampler& sampler){

#ifdef TEST_GLOBAL_TIMINGS
    kdTreeTime = 0;
    totalTime  = 0;
    verifyTime = 0;
#endif

    const Scalar kSmallError = 0.00001;
    const int kMinNumberOfTrials = 4;
    const Scalar kDiameterFraction = 0.3;

    centroid_P_ = VectorType::Zero();
    centroid_Q_ = VectorType::Zero();

	cv::Mat Pmat;
	//cv::Mat Pmat1;
	//P.pos();
    //Pmat = Mat(P[1].pos(), true);
	int pp;
	
	pp= P.size();
		const int DIM = pp;
	Mat indices(1, DIM, CV_32FC1);
	Mat dists(1, DIM, CV_32FC1);
	static float **kkk ;
	kkk = (float **)malloc(pp * sizeof(float *));
	for (int i = 0; i < pp; i++)
	{
		kkk[i] = (float *)malloc(3 * sizeof(float));
	}
	for (int i = 0; i < pp; i++)
	{
		kkk[i][0] = P[i].pos()[0];
		kkk[i][1] = P[i].pos()[1];
		kkk[i][2] = P[i].pos()[2];
	}

	//float KK[P.size()][3] = P[.pos();
	//linear_search(Pmat, Pmat, indices, dists, 8);
	//float m0[10][3] = { { 1,2,3},{6, 5, 4},	{7, 8, 9 },{10, 11, 12},{13, 14, 15},{16, 17, 18},{19, 20, 21},{22, 23, 24},{25, 26, 27},{28,29,30}};
	int m = sizeof(kkk[0]) / sizeof(float);
	int n = (sizeof(kkk) / sizeof(float)) / (sizeof(kkk[0]) / sizeof(float));
	//int mm = sizeof(kkkk[0]) / sizeof(int);
	//int nn = (sizeof(kkkk) / sizeof(int)) / (sizeof(kkkk[0]) / sizeof(int));
	Log<LogLevel::Verbose>("kkk长度: ", m,",kkk的宽度：",n);
	Log<LogLevel::Verbose>("kkk规格: ", sizeof(kkk));
	Log<LogLevel::Verbose>("pp: ", pp);

	Log<LogLevel::Verbose>("kkk[0,0]: ", kkk[0][0]);
	Log<LogLevel::Verbose>("P[i].pos()[0]: ", P[0].pos()[0]);
	Log<LogLevel::Verbose>("kkk[0,1]: ", kkk[0][1]);
	Log<LogLevel::Verbose>("P[0].pos()[1]: ", P[0].pos()[1]);
	Log<LogLevel::Verbose>("kkk[0,2]: ", kkk[0][2]);
	Log<LogLevel::Verbose>("P[0].pos()[2]: ", P[0].pos()[2]);
	Log<LogLevel::Verbose>("kkk[1,0]: ", kkk[1][0]);
	Log<LogLevel::Verbose>("P[1].pos()[0]: ", P[1].pos()[0]);
	Log<LogLevel::Verbose>("kkk[1,1]: ", kkk[1][1]);
	Log<LogLevel::Verbose>("P[1].pos()[1]: ", P[1].pos()[1]);
	Log<LogLevel::Verbose>("kkk[1,2]: ", kkk[1][2]);
	Log<LogLevel::Verbose>("P[1].pos()[2]: ", P[1].pos()[2]);
	Log<LogLevel::Verbose>("kkk[2,0]: ", kkk[2][0]);
	Log<LogLevel::Verbose>("P[2].pos()[0]: ", P[2].pos()[0]);
	Log<LogLevel::Verbose>("kkk[2,1]: ", kkk[2][1]);
	Log<LogLevel::Verbose>("P[2].pos()[1]: ", P[2].pos()[1]);
	Log<LogLevel::Verbose>("kkk[2,2]: ", kkk[2][2]);
	Log<LogLevel::Verbose>("P[2].pos()[2]: ", P[2].pos()[2]);
	//float aaa[100000][3];
	//static float aa[80000*3];
	//static float * aa = new float[80000 * 3];
	float * aa=new float[80000 * 3];
	for (int i = 0; i < pp; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			aa[(i-1)*3+j] = kkk[i][j];
		//	aaa[i][j] = aa;
		}
	}
	Log<LogLevel::Verbose>("aa[0]: ", aa[0]);
	Log<LogLevel::Verbose>("aa[1]: ", aa[1]);
	Log<LogLevel::Verbose>("aa[2]: ", aa[2]);
		cv::Mat mK = cv::Mat(5000,3, CV_32FC1, aa);

	    Log<LogLevel::Verbose>("P.size(): ", P.size());
	    Log<LogLevel::Verbose>("normals1.size(): ", normals1.size());
		Log<LogLevel::Verbose>("mK.size(): ", mK.size());

		clock_t ttt1 = clock();
		flann::Index flannIndex(mK, flann::KDTreeIndexParams());
		//flannIndex.buildIndex();
		//flann::Index flannIndex(mK, flann::LinearIndexParams());
		//flann::Index flannIndex(mK, flann::CompositeIndexParams());
		//flann::Index flannIndex(mK, flann::AutotunedIndexParams());
		clock_t ttt2 = clock();
		cout << "Begin LinearSearch: " << endl;

		clock_t tt1 = clock();
		flannIndex.knnSearch(mK, indices, dists, 8, flann::SearchParams());
		clock_t tt2 = clock();
		Log<LogLevel::Verbose>("索引时间: ", (ttt2 - ttt1) / CLOCKS_PER_SEC);
		Log<LogLevel::Verbose>("knn时间: ", (tt2 - tt1) / CLOCKS_PER_SEC);
		Log<LogLevel::Verbose>("knn+索引时间: ", (tt2 - ttt1) / CLOCKS_PER_SEC);

		Log<LogLevel::Verbose>("indices第一行: ", indices.row(1));
		Log<LogLevel::Verbose>("dists第一行: ", dists.row(1));
		Log<LogLevel::Verbose>("indices第二行: ", indices.row(2));
		Log<LogLevel::Verbose>("dists第二行: ", dists.row(2));
		
		//cout << "Finish linerSearch: " << (tt2 - tt1) / CLOCKS_PER_SEC << "s" << endl;
		
		//Log<LogLevel::Verbose>("point.size: ", point.size());
		Log<LogLevel::Verbose>("indices.size: ", indices.size());
		Log<LogLevel::Verbose>("dists.size: ", dists.size());
		//cout << indices << endl;
		
    sampled_P_3D_.clear();
    sampled_Q_3D_.clear();

    // prepare P
    if (P.size() > options_.sample_size){
        sampler(P, options_, sampled_P_3D_);
		Log<LogLevel::Verbose>("sampled_P_3D_size: ", sampled_P_3D_.size());
		//Pmat = Mat(sampled_P_3D_, true);
    }

    else
    {
        Log<LogLevel::ErrorReport>( "(P) More samples requested than available: use whole cloud" );
        sampled_P_3D_ = P;
    }



    // prepare Q
    if (Q.size() > options_.sample_size){
		clock_t k1 = clock();
        std::vector<Point3D> uniform_Q;
        sampler(Q, options_, uniform_Q);
		
	

        std::shuffle(uniform_Q.begin(), uniform_Q.end(), randomGenerator_);
        size_t nbSamples = std::min(uniform_Q.size(), options_.sample_size);
        auto endit = uniform_Q.begin(); std::advance(endit, nbSamples );
        std::copy(uniform_Q.begin(), endit, std::back_inserter(sampled_Q_3D_));
		clock_t k2 = clock();
		Log<LogLevel::Verbose>("uniform_Qsize: ", uniform_Q.size());
		Log<LogLevel::Verbose>("sampled_Q_3D_size: ", sampled_Q_3D_.size());
		Log<LogLevel::Verbose>("downsampled time: ", (k2-k1)/ CLOCKS_PER_SEC,"s");

    }
    else
    {
        Log<LogLevel::ErrorReport>( "(Q) More samples requested than available: use whole cloud" );
        sampled_Q_3D_ = Q;
    }


    // center points around centroids
    auto centerPoints = [](std::vector<Point3D>&container,
            VectorType& centroid){
        for(const auto& p : container) centroid += p.pos();
        centroid /= Scalar(container.size());
        for(auto& p : container) p.pos() -= centroid;
    };
    centerPoints(sampled_P_3D_, centroid_P_);
    centerPoints(sampled_Q_3D_, centroid_Q_);

    initKdTree();
    // Compute the diameter of P approximately (randomly). This is far from being
    // Guaranteed close to the diameter but gives good results for most common
    // objects if they are densely sampled.
    P_diameter_ = 0.0;
    for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
        int at = randomGenerator_() % sampled_Q_3D_.size();
        int bt = randomGenerator_() % sampled_Q_3D_.size();

        Scalar l = (sampled_Q_3D_[bt].pos() - sampled_Q_3D_[at].pos()).norm();
        if (l > P_diameter_) {
            P_diameter_ = l;
        }
    }

    // Mean distance and a bit more... We increase the estimation to allow for
    // noise, wrong estimation and non-uniform sampling.
    P_mean_distance_ = MeanDistance();

    // Normalize the delta (See the paper) and the maximum base distance.
    // delta = P_mean_distance_ * delta;
    max_base_diameter_ = P_diameter_;  // * estimated_overlap_;

    // RANSAC probability and number of needed trials.
    Scalar first_estimation =
            std::log(kSmallError) / std::log(1.0 - pow(options_.getOverlapEstimation(),
                                             static_cast<Scalar>(kMinNumberOfTrials)));
    // We use a simple heuristic to elevate the probability to a reasonable value
    // given that we don't simply sample from P, but instead, we bound the
    // distance between the points in the base as a fraction of the diameter.
    number_of_trials_ =
            static_cast<int>(first_estimation * (P_diameter_ / kDiameterFraction) /
                             max_base_diameter_);
    if (number_of_trials_ < kMinNumberOfTrials)
        number_of_trials_ = kMinNumberOfTrials;

    Log<LogLevel::Verbose>( "norm_max_dist: ", options_.delta );
    current_trial_ = 0;
    best_LCP_ = 0.0;

    Q_copy_ = Q;
    for (int i = 0; i < 4; ++i) {
        base_[i] = 0;
        current_congruent_[i] = 0;
    }
    transform_ = Eigen::Matrix<Scalar, 4, 4>::Identity();

    // call Virtual handler
    Initialize(P, Q);

    best_LCP_ = Verify(transform_);
    Log<LogLevel::Verbose>( "Initial LCP: ", best_LCP_ );
}


// Performs N RANSAC iterations and compute the best transformation. Also,
// transforms the set Q by this optimal transformation.
template <typename Visitor>
bool
Match4PCSBase::Perform_N_steps(int n,
                               Eigen::Ref<MatrixType> transformation,
                               std::vector<Point3D>* Q,
                               const Visitor &v) {
  using std::chrono::system_clock;
  if (Q == nullptr) return false;

#ifdef TEST_GLOBAL_TIMINGS
    Timer t (true);
#endif


  // The transformation has been computed between the two point clouds centered
  // at the origin, we need to recompute the translation to apply it to the original clouds
  auto getGlobalTransform = [this](Eigen::Ref<MatrixType> transformation){
    Eigen::Matrix<Scalar, 3, 3> rot, scale;
    Eigen::Transform<Scalar, 3, Eigen::Affine> (transform_).computeRotationScaling(&rot, &scale);
    transformation = transform_;
    transformation.col(3) = (qcentroid1_ + centroid_P_ - ( rot * scale * (qcentroid2_ + centroid_Q_))).homogeneous();
  };

  Scalar last_best_LCP = best_LCP_;
  v(0, best_LCP_, transformation);

  bool ok = false;
  clock_t start1, finish1;
  start1 = clock();
  //std::chrono::time_point<system_clock> t0 = system_clock::now(), end;
  for (int i = current_trial_; i < current_trial_ + n; ++i) {
    ok = TryOneBase(v);

    Scalar fraction_try  = Scalar(i) / Scalar(number_of_trials_);
	finish1 = clock();
    Scalar fraction_time =
      // std::chrono::duration_cast<std::chrono::seconds>
      //  (system_clock::now() - t0).count() /
		(finish1 - start1) / CLOCKS_PER_SEC/ options_.max_time_seconds;
    Scalar fraction = std::max(fraction_time, fraction_try);
	Log<LogLevel::Verbose>("循环: ", i);
	Log<LogLevel::Verbose>("用时: ", (finish1 - start1) / CLOCKS_PER_SEC);
    if (v.needsGlobalTransformation()) {
      getGlobalTransform(transformation);
    } else {
      transformation = transform_;
    }

    v(fraction, best_LCP_, transformation);

    // ok means that we already have the desired LCP.
	if (ok || i > number_of_trials_ || fraction >= 0.99 || best_LCP_ == 1.0)
	{
		//Log<LogLevel::Verbose>("ok: ", ok);
		//Log<LogLevel::Verbose>("i: ", i);
		//Log<LogLevel::Verbose>("number_of_trials_: ", number_of_trials_);
		//Log<LogLevel::Verbose>("fraction_time: ", fraction_time);
		//Log<LogLevel::Verbose>("fraction_try: ", fraction_try);
		//Log<LogLevel::Verbose>("fraction: ", fraction);
		break;
	}
  }

  current_trial_ += n;
  if (best_LCP_ > last_best_LCP) {
    *Q = Q_copy_;

    getGlobalTransform(transformation);

    // Transforms Q by the new transformation.
    for (size_t i = 0; i < Q->size(); ++i) {
      (*Q)[i].pos() = (transformation * (*Q)[i].pos().homogeneous()).head<3>();
    }
  }
#ifdef TEST_GLOBAL_TIMINGS
    totalTime += Scalar(t.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif

  return ok || current_trial_ >= number_of_trials_;
}



// Pick one base, finds congruent 4-points in Q, verifies for all
// transformations, and retains the best transformation and LCP. This is
// a complete RANSAC iteration.
template<typename Visitor>
bool Match4PCSBase::TryOneBase(const Visitor &v) {
  Scalar invariant1, invariant2;
  int base_id1, base_id2, base_id3, base_id4;

//#define STATIC_BASE

#ifdef STATIC_BASE
  static bool first_time = true;

  if (first_time){
      base_id1 = 0;
      base_id2 = 3;
      base_id3 = 1;
      base_id4 = 4;

      base_3D_[0] = sampled_P_3D_ [base_id1];
      base_3D_[1] = sampled_P_3D_ [base_id2];
      base_3D_[2] = sampled_P_3D_ [base_id3];
      base_3D_[3] = sampled_P_3D_ [base_id4];

      TryQuadrilateral(&invariant1, &invariant2, base_id1, base_id2, base_id3, base_id4);

      first_time = false;
  }
  else
      return false;

#else

  if (!SelectQuadrilateral(invariant1, invariant2, base_id1, base_id2,
                           base_id3, base_id4)) {
    return false;
  }
#endif

  // Computes distance between pairs.
  const Scalar distance1 = (base_3D_[0].pos()- base_3D_[1].pos()).norm();
  const Scalar distance2 = (base_3D_[2].pos()- base_3D_[3].pos()).norm();

  std::vector<std::pair<int, int>> pairs1, pairs2;
  std::vector<Quadrilateral> congruent_quads;

  // Compute normal angles.
  const Scalar normal_angle1 = (base_3D_[0].normal() - base_3D_[1].normal()).norm();
  const Scalar normal_angle2 = (base_3D_[2].normal() - base_3D_[3].normal()).norm();

  ExtractPairs(distance1, normal_angle1, distance_factor * options_.delta, 0,
                  1, &pairs1);
  ExtractPairs(distance2, normal_angle2, distance_factor * options_.delta, 2,
                  3, &pairs2);

//  Log<LogLevel::Verbose>( "Pair creation ouput: ", pairs1.size(), " - ", pairs2.size());

  if (pairs1.size() == 0 || pairs2.size() == 0) {
    return false;
  }


  if (!FindCongruentQuadrilaterals(invariant1, invariant2,
                                   distance_factor * options_.delta,
                                   distance_factor * options_.delta,
                                   pairs1,
                                   pairs2,
                                   &congruent_quads)) {
    return false;
  }

  size_t nb = 0;

  bool match = TryCongruentSet(base_id1, base_id2, base_id3, base_id4,
                               congruent_quads,
                               v,
                               nb);

  //if (nb != 0)
  //  Log<LogLevel::Verbose>( "Congruent quads: (", nb, ")    " );

  return match;
}


template <typename Visitor>
bool Match4PCSBase::TryCongruentSet(
        int base_id1,
        int base_id2,
        int base_id3,
        int base_id4,
        const std::vector<Quadrilateral>& congruent_quads,
        const Visitor& v,
        size_t &nbCongruent){
    static const double pi = std::acos(-1);

    // get references to the basis coordinates
    const Point3D& b1 = sampled_P_3D_[base_id1];
    const Point3D& b2 = sampled_P_3D_[base_id2];
    const Point3D& b3 = sampled_P_3D_[base_id3];
    const Point3D& b4 = sampled_P_3D_[base_id4];

    // set the basis coordinates in the congruent quad array
    const std::array<Point3D, 4> congruent_base {{b1, b2, b3, b4}};


    // Centroid of the basis, computed once and using only the three first points
    Eigen::Matrix<Scalar, 3, 1> centroid1 = (b1.pos() + b2.pos() + b3.pos()) / Scalar(3);


    std::atomic<size_t> nbCongruentAto(0);

#ifdef SUPER4PCS_USE_OPENMP
#pragma omp parallel for num_threads(omp_nthread_congruent_)
#endif
    for (int i = 0; i < int(congruent_quads.size()); ++i) {
      std::array<Point3D, 4> congruent_candidate;

      Eigen::Matrix<Scalar, 4, 4> transform;

      // Centroid of the sets, computed in the loop using only the three first points
      Eigen::Matrix<Scalar, 3, 1> centroid2;

      const int a = congruent_quads[i].vertices[0];
      const int b = congruent_quads[i].vertices[1];
      const int c = congruent_quads[i].vertices[2];
      const int d = congruent_quads[i].vertices[3];
      congruent_candidate[0] = sampled_Q_3D_[a];
      congruent_candidate[1] = sampled_Q_3D_[b];
      congruent_candidate[2] = sampled_Q_3D_[c];
      congruent_candidate[3] = sampled_Q_3D_[d];

  #ifdef STATIC_BASE
      Log<LogLevel::Verbose>( "Ids: ", base_id1, "\t", base_id2, "\t", base_id3, "\t", base_id4);
      Log<LogLevel::Verbose>( "     ", a, "\t", b, "\t", c, "\t", d);
  #endif

      centroid2 = (congruent_candidate[0].pos() +
                   congruent_candidate[1].pos() +
                   congruent_candidate[2].pos()) / Scalar(3.);

      Scalar rms = -1;

      const bool ok =
      ComputeRigidTransformation(congruent_base,     // input congruent quad
                                 congruent_candidate,// tested congruent quad
                                 centroid1,          // input: basis centroid
                                 centroid2,          // input: candidate quad centroid
                                 options_.max_angle * pi / 180.0, // maximum per-dimension angle, check return value to detect invalid cases
                                 transform,          // output: transformation
                                 rms,                // output: rms error of the transformation between the basis and the congruent quad
                             #ifdef MULTISCALE
                                 true
                             #else
                                 false
                             #endif
                                 );             // state: compute scale ratio ?

      if (ok && rms >= Scalar(0.)) {

        // We give more tolerantz in computing the best rigid transformation.
        if (rms < distance_factor * options_.delta) {

          nbCongruentAto++;
          // The transformation is computed from the point-clouds centered inn [0,0,0]

          // Verify the rest of the points in Q against P.
          Scalar lcp = Verify(transform);

          // transformation has been computed between the two point clouds centered
          // at the origin, we need to recompute the translation to apply it to the original clouds
          auto getGlobalTransform =
              [this, transform, centroid1, centroid2]
              (Eigen::Ref<MatrixType> transformation){
            Eigen::Matrix<Scalar, 3, 3> rot, scale;
            Eigen::Transform<Scalar, 3, Eigen::Affine> (transform).computeRotationScaling(&rot, &scale);
            transformation = transform;
            transformation.col(3) = (centroid1 + centroid_P_ - ( rot * scale * (centroid2 + centroid_Q_))).homogeneous();
          };

          if (v.needsGlobalTransformation())
          {
            Eigen::Matrix<Scalar, 4, 4> transformation = transform;
            getGlobalTransform(transformation);
            v(-1, lcp, transformation);
          }
          else
            v(-1, lcp, transform);

#pragma omp critical
          if (lcp > best_LCP_) {
            // Retain the best LCP and transformation.
            base_[0] = base_id1;
            base_[1] = base_id2;
            base_[2] = base_id3;
            base_[3] = base_id4;

            current_congruent_[0] = a;
            current_congruent_[1] = b;
            current_congruent_[2] = c;
            current_congruent_[3] = d;

            best_LCP_    = lcp;
            transform_   = transform;
            qcentroid1_  = centroid1;
            qcentroid2_  = centroid2;
          }
          // Terminate if we have the desired LCP already.
          if (best_LCP_ > options_.getTerminateThreshold()){
            continue;
          }
        }
      }
    }

    nbCongruent = nbCongruentAto;

    // If we reached here we do not have yet the desired LCP.
    return best_LCP_ > options_.getTerminateThreshold() /*false*/;
}


} // namespace Super4PCS

#endif
