#include "super4pcs/algorithms/4pcs.h"
#include "super4pcs/algorithms/super4pcs.h"
#include "super4pcs/io/io.h"
#include "super4pcs/utils/geometry.h"
#include<ctime>
#include <E:\Program Files\PCL 1.8.1\3rdParty\Eigen\eigen3\Eigen/Dense>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "../demo-utils.h"
#define sqr(x) ((x) * (x))

using namespace std;
using namespace GlobalRegistration;
using namespace GlobalRegistration::Demo;


static inline void printS4PCSParameterList(){
	fprintf(stderr, "\t[ Welcome to use Super4pcs revised by SIA Rancho ]\n");
    fprintf(stderr, "\t[ -r result_file_name (%s) ]\n", output.c_str());
    fprintf(stderr, "\t[ -m output matrix file (%s) ]\n", outputMat.c_str());
    fprintf(stderr, "\t[ -x (use 4pcs: false by default) ]\n");
    fprintf(stderr, "\t[ --sampled1 (output sampled cloud 1 -- debug+super4pcs only) ]\n");
    fprintf(stderr, "\t[ --sampled2 (output sampled cloud 2 -- debug+super4pcs only) ]\n");
}

struct TransformVisitor {
    inline void operator()(
            float fraction,
            float best_LCP,
            Eigen::Ref<Match4PCSBase::MatrixType> /*transformation*/) const {
      if (fraction >= 0)
        {
          printf("done: %d%c      best: %f                  \r",
               static_cast<int>(fraction * 100), '%', best_LCP);
          fflush(stdout);
        }
    }
    constexpr bool needsGlobalTransformation() const { return false; }
};
double _f(double t, int n) {
	double m = pow(10, n);
	double result = floor(t * m + 0.5) / m;
	return result;
}
float youxiaoshuzi(float num)
{
	double result;
	float num1 = num;
	int q = 0;
	while(1)
	{
	if (floor(num1) > 1)
	{
	result=_f(num, q);
		break;
	}
	q = q + 1;
	num1 = num1 * 10; 
	}
	return result;
}
int main(int argc, char **argv) {
  using namespace GlobalRegistration;

  vector<Point3D> set1, set2;
  vector<Eigen::Matrix2f> tex_coords1, tex_coords2;
  vector<typename Point3D::VectorType> normals1, normals2;
  vector<tripple> tris1, tris2;
  vector<std::string> mtls1, mtls2;

  // Match and return the score (estimated overlap or the LCP).
  typename Point3D::Scalar score = 0;

  constexpr Utils::LogLevel loglvl = Utils::Verbose;
  using SamplerType   = GlobalRegistration::Sampling::UniformDistSampler;
  using TrVisitorType = typename std::conditional <loglvl==Utils::NoLog,
                            Match4PCSBase::DummyTransformVisitor,
                            TransformVisitor>::type;
  SamplerType sampler;
  TrVisitorType visitor;
  Utils::Logger logger(loglvl);

  /// TODO Add proper error codes
  if(argc < 4){
      Demo::printUsage(argc, argv);
      exit(-2);
  }
  if(int c = Demo::getArgs(argc, argv) != 0)
  {
    Demo::printUsage(argc, argv);
    printS4PCSParameterList();
    exit(std::max(c,0));
  }

  Match4PCSOptions options;

  Match4PCSBase::MatrixType mat (Match4PCSBase::MatrixType::Identity());

  if(! Demo::setOptionsFromArgs(options, logger))   //这里把options的数据重新刷新了一遍
  {
    exit(-3);
  }
  //logger.Log<Utils::Verbose>("options.delta3:", options.delta);
  // load data
  IOManager iomananger;
  clock_t start, finish;
  start = clock();
  // Read the inputs.
  if (!iomananger.ReadObject((char *)input1.c_str(), set1, tex_coords1, normals1, tris1,
                  mtls1)) {
    logger.Log<Utils::ErrorReport>("Can't read input set1");
    exit(-1);
  }
  int aa = set1.size();
  float qxmax= set1[0].pos()[0];
  float qxmin= set1[0].pos()[0];
  float qymax = set1[0].pos()[1];
  float qymin = set1[0].pos()[1];
  float qzmax = set1[0].pos()[2];
  float qzmin = set1[0].pos()[2];
  for (int i=0;i<aa-1;i++)
  {
	  float qx = set1[i].pos()[0];
	  qxmax = max(qx,qxmax);
	  qxmin = min(qx,qxmin);
	  float qy = set1[i].pos()[1];
	  qymax = max(qy,qymax);
	  qymin = min(qy,qymin);
	  float qz = set1[i].pos()[2];
	  qzmax = max(qz,qzmax);
	  qzmin = min(qz,qzmin);

}  
  //logger.Log<Utils::Verbose>("点云1的法向量:", normals1[1]);
  //logger.Log<Utils::Verbose>("点云1的法向量:", normals1[1][0]);
  //logger.Log<Utils::Verbose>("点云1的法向量:", normals1[1][1]);
  //logger.Log<Utils::Verbose>("点云1的法向量:", normals1[1][2]);
  //logger.Log<Utils::Verbose>("点云1的x坐标最大值:", qxmax);
  //logger.Log<Utils::Verbose>("点云1的x坐标最小值:", qxmin);
  //logger.Log<Utils::Verbose>("点云1的y坐标最大值:", qymax);
  //logger.Log<Utils::Verbose>("点云1的y坐标最小值:", qymin);
  //logger.Log<Utils::Verbose>("点云1的z坐标最大值:", qzmax);
  //logger.Log<Utils::Verbose>("点云1的z坐标最小值:", qzmin);
  double pnum = (set1.size());
  logger.Log<Utils::Verbose>("点云1的点数量:", pnum);
  float pointnumber = pow(pnum,1.0/3);
  logger.Log<Utils::Verbose>("pointnumber:", pointnumber);
  //logger.Log<Utils::Verbose>("点云1的x坐标值范围:", (qxmax - qxmin));
  //logger.Log<Utils::Verbose>("点云1的y坐标值范围:", (qymax - qymin));
  //logger.Log<Utils::Verbose>("点云1的z坐标值范围:", (qzmax - qzmin));
  logger.Log<Utils::Verbose>("时间限制:", max_time_seconds);
  float deltamax = max(max((qxmax - qxmin) / pointnumber, (qymax - qymin) / pointnumber), (qzmax - qzmin) / pointnumber);
  logger.Log<Utils::Verbose>("计算的delta:", deltamax);
  logger.Log<Utils::Verbose>("options.sample_size:", options.sample_size);
  delta = youxiaoshuzi(deltamax);
  options.delta = youxiaoshuzi(deltamax);
  logger.Log<Utils::Verbose>("设置的delta:", delta);
  logger.Log<Utils::Verbose>("设置的options.delta:", options.delta);
  if (!iomananger.ReadObject((char *)input2.c_str(), set2, tex_coords2, normals2, tris2,
                  mtls2)) {
    logger.Log<Utils::ErrorReport>("Can't read input set2");
    exit(-1);
  }
  logger.Log<Utils::Verbose>("点云2的点数量:", set2.size());
  // clean only when we have pset to avoid wrong face to point indexation
  if (tris1.size() == 0)
    Utils::CleanInvalidNormals(set1, normals1);
  if (tris2.size() == 0)
    Utils::CleanInvalidNormals(set2, normals2);

  try {

      if (use_super4pcs) {
          MatchSuper4PCS matcher(options, logger);  //这句话定义在哪
          logger.Log<Utils::Verbose>( "Use Super4PCS revised by SIA Rancho" );
	          score = matcher.ComputeTransformation(set1,normals1, &set2,&normals2, mat, sampler, visitor);

          if(! outputSampled1.empty() ){
              logger.Log<Utils::Verbose>( "Exporting Sampled cloud 1 to ",
                                          outputSampled1.c_str(),
                                          " ..." );
              iomananger.WriteObject((char *)outputSampled1.c_str(),
                                     matcher.getFirstSampled(),
                                     vector<Eigen::Matrix2f>(),
                                     vector<typename Point3D::VectorType>(),
                                     vector<tripple>(),
                                     vector<string>());
              logger.Log<Utils::Verbose>( "Export DONE" );
          }
          if(! outputSampled2.empty() ){
              logger.Log<Utils::Verbose>( "Exporting Sampled cloud 2 to ",
                                          outputSampled2.c_str(),
                                          " ..." );
              iomananger.WriteObject((char *)outputSampled2.c_str(),
                                     matcher.getSecondSampled(),
                                     vector<Eigen::Matrix2f>(),
                                     vector<typename Point3D::VectorType>(),
                                     vector<tripple>(),
                                     vector<string>());
              logger.Log<Utils::Verbose>( "Export DONE" );
          }
      }
      else {
          Match4PCS matcher(options, logger);
          logger.Log<Utils::Verbose>( "Use old 4PCS" );
          score = matcher.ComputeTransformation(set1, normals1, &set2, &normals2, mat, sampler, visitor);
      }

  }
  catch (const std::exception& e) {
      logger.Log<Utils::ErrorReport>( "[Error]: " , e.what() );
      logger.Log<Utils::ErrorReport>( "Aborting with code -2 ..." );
      return -2;
  }
  catch (...) {
      logger.Log<Utils::ErrorReport>( "[Unknown Error]: Aborting with code -3 ..." );
      return -3;
  }
  finish = clock();
  logger.Log<Utils::Verbose>("Total Time: ", (finish - start)/ CLOCKS_PER_SEC,"s");
  logger.Log<Utils::Verbose>( "Score: ", score );
  logger.Log<Utils::Verbose>( "(Homogeneous) Transformation from ",
                              input2.c_str(),
                              " to ",
                              input1.c_str(),
                              ": \n",
                              mat);


  if(! outputMat.empty() ){
      logger.Log<Utils::Verbose>( "Exporting Matrix to ",
                                  outputMat.c_str(),
                                  "..." );
      iomananger.WriteMatrix(outputMat, mat.cast<double>(), IOManager::POLYWORKS);
      logger.Log<Utils::Verbose>( "Export DONE" );
  }

  if (! output.empty() ){

      logger.Log<Utils::Verbose>( "Exporting Registered geometry to ",
                                  output.c_str(),
                                  "..." );
      iomananger.WriteObject((char *)output.c_str(),
                             set2,
                             tex_coords2,
                             normals2,
                             tris2,
                             mtls2);
      logger.Log<Utils::Verbose>( "Export DONE" );
  }

  return 0;
}
