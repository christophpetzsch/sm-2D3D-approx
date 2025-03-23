#ifndef ARAP_LOSS
#define ARAP_LOSS

#include <Eigen/Dense>

#define DEBUG_ARAP false    
#define EPSI 1e-8

#define AVG_NORMALS true
#define PRECOMPUTE_ROTATIONS true
#define EPSI 1e-8

const double getLocalThickness3d(const int vertexIdx, const double* features3d, const int dimFeatures3d) {
    return features3d[dimFeatures3d * vertexIdx + 14];
}

const double getLocalThickness2d(const int vertexIdx, const double* features2d, const int dimFeatures2d) {
    return features2d[dimFeatures2d * vertexIdx + 7];
}

inline double robustLoss(const double dist2Other2d, const double dist2Other3d) {
    const double alpha = -2;
    const double c = 0.15; 
    const double x = dist2Other2d - dist2Other3d;
    
    const double xDivCSqr = std::pow((x / c), 2);
    const double absAlphaMinusTwo = std::abs(alpha - 2);
    
    return absAlphaMinusTwo / alpha * (std::pow(xDivCSqr / absAlphaMinusTwo + 1, 0.5 * alpha) -1);
}

double arapLoss(const double*   features3d,
                const int       nVertices3d,
                const int       dimFeatures3d,
                const double*   features2d,
                const int       nVertices2d,
                const int       dimFeatures2d,
                const long      current2dSrc,
                const long      current2dTrgt,
                const long      next2dSrc,
                const long      next2dTrgt,
                const long      current3dSrc,
                const long      current3dTrgt,
                const long      next3dSrc,
                const long      next3dTrgt) {
    if (DEBUG_ARAP) {
        if (current2dSrc < 0 || current2dSrc > nVertices2d)
            std::cout << "current2dSrc out of bounds" << std::endl;
        if (current2dTrgt < 0 || current2dTrgt > nVertices2d)
            std::cout << "current2dTrgt out of bounds" << std::endl;
        if (next2dSrc < 0 || next2dSrc > nVertices2d)
            std::cout << "next2dSrc out of bounds" << std::endl;
        if (next2dTrgt < 0 || next2dTrgt > nVertices2d)
            std::cout << "next2dTrgt out of bounds" << std::endl;
        
        if (current2dSrc == current2dTrgt)
            std::cout << "current2dSrc == current2dTrgt" << std::endl;
        if (next2dSrc == next2dTrgt)
            std::cout << "next2dSrc == next2dTrgt" << std::endl;
        if (current2dSrc+1 != current2dTrgt)
            std::cout << "current2dSrc("<< current2dSrc <<")+1 != current2dTrgt("<<current2dTrgt<<") " << std::endl;
        if (next2dSrc+1 != next2dTrgt)
            std::cout << "next2dSrc("<< next2dSrc <<")+1!= next2dTrgt("<<next2dTrgt<<") " << std::endl;
    }
    
    
    const Eigen::Vector3d v1_2d(features2d[current2dTrgt*dimFeatures2d  ] - features2d[current2dSrc*dimFeatures2d],
                                features2d[current2dTrgt*dimFeatures2d+1] - features2d[current2dSrc*dimFeatures2d+1],
                                features2d[current2dTrgt*dimFeatures2d+2] - features2d[current2dSrc*dimFeatures2d+2]);
    const Eigen::Vector3d v2_2d(features2d[next2dTrgt*dimFeatures2d  ] - features2d[next2dSrc*dimFeatures2d],
                                features2d[next2dTrgt*dimFeatures2d+1] - features2d[next2dSrc*dimFeatures2d+1],
                                features2d[next2dTrgt*dimFeatures2d+2] - features2d[next2dSrc*dimFeatures2d+2]);
    
    if (DEBUG_ARAP && v1_2d.norm() < EPSI) std::cout << "v1_2d zero" << std::endl;
    if (DEBUG_ARAP && v2_2d.norm() < EPSI) std::cout << "v2_2d zero" << std::endl;
    
    const Eigen::Vector3d v1_3d(features3d[current3dTrgt*dimFeatures3d  ] - features3d[current3dSrc*dimFeatures3d],
                                features3d[current3dTrgt*dimFeatures3d+1] - features3d[current3dSrc*dimFeatures3d+1],
                                features3d[current3dTrgt*dimFeatures3d+2] - features3d[current3dSrc*dimFeatures3d+2]);
    const Eigen::Vector3d v2_3d(features3d[next3dTrgt*dimFeatures3d  ] - features3d[next3dSrc*dimFeatures3d],
                                features3d[next3dTrgt*dimFeatures3d+1] - features3d[next3dSrc*dimFeatures3d+1],
                                features3d[next3dTrgt*dimFeatures3d+2] - features3d[next3dSrc*dimFeatures3d+2]);
    
    if (DEBUG_ARAP && v1_3d.norm() < EPSI) std::cout << "v1_3d zero" << std::endl;
    if (DEBUG_ARAP && v2_3d.norm() < EPSI) std::cout << "v2_3d zero" << std::endl;
    
    const double dist2Other1_2d = getLocalThickness2d(current2dTrgt, features2d, dimFeatures2d);
    const double dist2Other1_3d = getLocalThickness3d(current3dTrgt, features3d, dimFeatures3d);
    const double dist2Other2_2d = getLocalThickness2d(next2dTrgt, features2d, dimFeatures2d);
    const double dist2Other2_3d = getLocalThickness3d(next3dTrgt, features3d, dimFeatures3d);
    
    
    Eigen::Vector3d n1_2d(features2d[current2dSrc * dimFeatures2d + 3], 
                                features2d[current2dSrc * dimFeatures2d + 4], 
                                features2d[current2dSrc * dimFeatures2d + 5]);
    if (DEBUG_ARAP && n1_2d.norm() < EPSI) std::cout << "n1_2d zero" << std::endl;
    Eigen::Vector3d n2_2d(features2d[next2dSrc * dimFeatures2d + 3], 
                                features2d[next2dSrc * dimFeatures2d + 4], 
                                features2d[next2dTrgt * dimFeatures2d + 5]);
    
    if (DEBUG_ARAP && n2_2d.norm() < EPSI) std::cout << "n2_2d zero" << std::endl;
    
    Eigen::Vector3d n1_3d(features3d[current3dSrc * dimFeatures3d + 11], 
                                features3d[current3dSrc * dimFeatures3d + 12], 
                                features3d[current3dSrc * dimFeatures3d + 13]);
    if (AVG_NORMALS) {
        Eigen::Vector3d n1_3dTemp(features3d[current3dTrgt * dimFeatures3d + 11], 
                                features3d[current3dTrgt * dimFeatures3d + 12], 
                                features3d[current3dTrgt * dimFeatures3d + 13]);
        n1_3d = (n1_3d + n1_3dTemp) / 2;
    }
    if (DEBUG_ARAP && n1_3d.norm() < EPSI) std::cout << "n1_3d zero" << std::endl;
    
    Eigen::Vector3d n2_3d(features3d[next3dSrc * dimFeatures3d + 11], 
                                features3d[next3dSrc * dimFeatures3d + 12], 
                                features3d[next3dSrc * dimFeatures3d + 13]);
    if (AVG_NORMALS) {
        Eigen::Vector3d n2_3dTemp(features3d[next3dTrgt * dimFeatures3d + 11], 
                                features3d[next3dTrgt * dimFeatures3d + 12], 
                                features3d[next3dTrgt * dimFeatures3d + 13]);
        n2_3d = (n2_3d + n2_3dTemp) / 2;
    }
    if (DEBUG_ARAP && n2_3d.norm() < EPSI) std::cout << "n2_3d zero" << std::endl;
    //n2_3d *= dist2Other2_3d;
    
    const Eigen::Vector3d c1_2d = v1_2d.cross(n1_2d);
    const Eigen::Vector3d c2_2d = v2_2d.cross(n2_2d);
    
    const Eigen::Vector3d c1_3d = v1_3d.cross(n1_3d);
    const Eigen::Vector3d c2_3d = v2_3d.cross(n2_3d);

    

    Eigen::Matrix3d X1;
    X1.row(0) = n1_2d.normalized();
    X1.row(1) = v1_2d.normalized();
    X1.row(2) = c1_2d.normalized();
    Eigen::Matrix3d X2;
    X2.row(0) = n2_2d.normalized();
    X2.row(1) = v2_2d.normalized();
    X2.row(2) = c2_2d.normalized();
    
    Eigen::Matrix3d Y1;
    Y1.row(0) = n1_3d.normalized();
    Y1.row(1) = v1_3d.normalized();
    Y1.row(2) = c1_3d.normalized();
    Eigen::Matrix3d Y2;
    Y2.row(0) = n2_3d.normalized();
    Y2.row(1) = v2_3d.normalized();
    Y2.row(2) = c2_3d.normalized();
    
    
    Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeFullV | Eigen::ComputeFullU > svd1( X1.transpose() * Y1);
    Eigen::DiagonalMatrix<double, 3> d1(1, 1, (svd1.matrixU() * (svd1.matrixV().transpose())).determinant() );


    Eigen::Matrix3d rot1 = svd1.matrixU() * d1 * (svd1.matrixV().transpose());
    
    Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeFullV | Eigen::ComputeFullU > svd2( X2.transpose() * Y2 );
    Eigen::DiagonalMatrix<double, 3> d2(1, 1, (svd2.matrixU() * (svd2.matrixV().transpose())).determinant() );
    
    if (DEBUG_ARAP && (svd1.singularValues().transpose().array().abs() < EPSI).any()) {
        std::cout << "svd1 singular values zero" << std::endl;
    }
    if (DEBUG_ARAP && (svd2.singularValues().transpose().array().abs() < EPSI).any()) {
        std::cout << "svd2 singular values zero" << std::endl;
    }
    
    Eigen::Matrix3d rot2 = svd2.matrixU() * d2 * (svd2.matrixV().transpose());
    
    if (DEBUG_ARAP && (X1.transpose() * Y1 - svd1.matrixU() * svd1.singularValues().asDiagonal() * svd1.matrixV().transpose()).norm() > EPSI)
            std::cout << "SVD1: A != USV^T" << std::endl;
    if (DEBUG_ARAP && (X2.transpose() * Y2 - svd2.matrixU() * svd2.singularValues().asDiagonal() * svd2.matrixV().transpose()).norm() > EPSI)
            std::cout << "SVD2: A != USV^T" << std::endl;
    
    const Eigen::Quaterniond rotation1(rot1);
    //std::cout << rotation1.w() << std::endl;
	const Eigen::Quaterniond rotation2(rot2);
    //std::cout << rotation2.w() << std::endl;
    
    double innerProductRot2 = fabs(rotation1.w()*rotation2.w() + rotation1.vec().dot(rotation2.vec()));
    if ( innerProductRot2 > 1 ) {
        innerProductRot2 = 1;	
    } 
    if ( innerProductRot2 < -1 ) { 
        innerProductRot2 = -1;
    }

    
    double angleDifference2 = 2*acos(innerProductRot2);
    return std::abs(angleDifference2);
}
/*
 *
 *
 *
 *
 */
Eigen::Quaterniond precomputeRotation(const double*   features3d,
                                        const int       nVertices3d,
                                        const int       dimFeatures3d,
                                        const double*   features2d,
                                        const int       nVertices2d,
                                        const int       dimFeatures2d,
                                        const long      current2dSrc,
                                        const long      current2dTrgt,
                                        const long      current3dSrc,
                                        const long      current3dTrgt) {
    if (DEBUG_ARAP) {
        if (current2dSrc < 0 || current2dSrc > nVertices2d)
            std::cout << "current2dSrc out of bounds" << std::endl;
        if (current2dTrgt < 0 || current2dTrgt > nVertices2d)
            std::cout << "current2dTrgt out of bounds" << std::endl;
        if (current2dSrc == current2dTrgt)
            std::cout << "current2dSrc == current2dTrgt" << std::endl;
        if (current2dSrc+1 != current2dTrgt)
            std::cout << "current2dSrc("<< current2dSrc <<")+1 != current2dTrgt("<<current2dTrgt<<") " << std::endl;
    }
    
    
    const Eigen::Vector3d v1_2d(features2d[current2dTrgt*dimFeatures2d  ] - features2d[current2dSrc*dimFeatures2d],
                                features2d[current2dTrgt*dimFeatures2d+1] - features2d[current2dSrc*dimFeatures2d+1],
                                features2d[current2dTrgt*dimFeatures2d+2] - features2d[current2dSrc*dimFeatures2d+2]);
    if (DEBUG_ARAP && v1_2d.norm() < EPSI) std::cout << "v1_2d zero" << std::endl;
    
    const Eigen::Vector3d v1_3d(features3d[current3dTrgt*dimFeatures3d  ] - features3d[current3dSrc*dimFeatures3d],
                                features3d[current3dTrgt*dimFeatures3d+1] - features3d[current3dSrc*dimFeatures3d+1],
                                features3d[current3dTrgt*dimFeatures3d+2] - features3d[current3dSrc*dimFeatures3d+2]);
    if (DEBUG_ARAP && v1_3d.norm() < EPSI) std::cout << "v1_3d zero" << std::endl;

    
    Eigen::Vector3d n1_2d(features2d[current2dSrc * dimFeatures2d + 3], 
                                features2d[current2dSrc * dimFeatures2d + 4], 
                                features2d[current2dSrc * dimFeatures2d + 5]);
    if (DEBUG_ARAP && n1_2d.norm() < EPSI) std::cout << "n1_2d zero" << std::endl;

    
    Eigen::Vector3d n1_3d(features3d[current3dSrc * dimFeatures3d + 11], 
                                features3d[current3dSrc * dimFeatures3d + 12], 
                                features3d[current3dSrc * dimFeatures3d + 13]);
    if (AVG_NORMALS) {
        Eigen::Vector3d n1_3dTemp(features3d[current3dTrgt * dimFeatures3d + 11], 
                                features3d[current3dTrgt * dimFeatures3d + 12], 
                                features3d[current3dTrgt * dimFeatures3d + 13]);
        n1_3d = (n1_3d + n1_3dTemp) / 2;
    }
    if (DEBUG_ARAP && n1_3d.norm() < EPSI) std::cout << "n1_3d zero" << std::endl;

    
    const Eigen::Vector3d c1_2d = v1_2d.cross(n1_2d);
    const Eigen::Vector3d c1_3d = v1_3d.cross(n1_3d);

    
    Eigen::Matrix3d X1;
    X1.row(0) = n1_2d.normalized();
    X1.row(1) = v1_2d.normalized();
    X1.row(2) = c1_2d.normalized();
 
    Eigen::Matrix3d Y1;
    Y1.row(0) = n1_3d.normalized();
    Y1.row(1) = v1_3d.normalized();
    Y1.row(2) = c1_3d.normalized();


    Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeFullV | Eigen::ComputeFullU > svd1( X1.transpose() * Y1);
    Eigen::DiagonalMatrix<double, 3> d1(1, 1, (svd1.matrixU() * (svd1.matrixV().transpose())).determinant() );


    Eigen::Matrix3d rot1 = svd1.matrixU() * d1 * (svd1.matrixV().transpose());

    if (DEBUG_ARAP && (svd1.singularValues().transpose().array().abs() < EPSI).any()) {
        std::cout << "svd1 singular values zero" << std::endl;
    }

    if (DEBUG_ARAP && (X1.transpose() * Y1 - svd1.matrixU() * svd1.singularValues().asDiagonal() * svd1.matrixV().transpose()).norm() > EPSI)
            std::cout << "SVD1: A != USV^T" << std::endl;
    
    return Eigen::Quaterniond(rot1);
}

/*
 *
 *
 *
 *
 */

inline double robustArapDif(const double arap) {
    const double alpha = -1.2;
    const double c = 0.2;

    const double x = arap;
    
    const double xDivCSqr = std::pow((x / c), 2);
    const double absAlphaMinusTwo = std::abs(alpha - 2);
    
    return absAlphaMinusTwo / alpha * (std::pow(xDivCSqr / absAlphaMinusTwo + 1, 0.5 * alpha) -1);    
}
/*
 *
 *
 *
 *
 */
double getArapWPrecompute(const long      numProductNodesPerLayer,
                          const long      currentIdxInLayer,
                          const long      nextIdxInLayer,
                          const long      current2dSrc,
                          const long      next2dSrc,
                          const std::vector< Eigen::Quaterniond>& precRotations) {
    const long idx1 = current2dSrc * numProductNodesPerLayer + currentIdxInLayer;
    const long idx2 = next2dSrc * numProductNodesPerLayer + nextIdxInLayer;
    const Eigen::Quaterniond rotation1 = precRotations.at(idx1);
	const Eigen::Quaterniond rotation2 = precRotations.at(idx2);
    double innerProductRot2 = fabs(rotation1.w()*rotation2.w() + rotation1.vec().dot(rotation2.vec()));
    if ( innerProductRot2 > 1 ) {
        innerProductRot2 = 1;	
    } 
    if ( innerProductRot2 < -1 ) { 
        innerProductRot2 = -1;
    }
    return 2 * acos(innerProductRot2);
}
/*
 *
 *
 *
 *
 */
inline double robustLossArap(const double x) {
    const double alpha = 0.7;
    const double c = 0.6;
    
    const double xDivCSqr = std::pow((x / c), 4);
    const double absAlphaMinusTwo = std::abs(alpha - 2);
    
    return absAlphaMinusTwo / alpha * (std::pow(xDivCSqr / absAlphaMinusTwo + 1, 0.5 * alpha) -1);    
}
/*
 *
 *
 *
 *
 */
inline double getArapCost (        
        const double*   features3d,
        const int       nVertices3d,
        const int       dimFeatures3d,
        const double*   features2d,
        const int       nVertices2d,
        const int       dimFeatures2d,
        const long      currentIdxInLayer,
        const long      nextIdxInLayer,
        const long      current2dSrc,
        const long      next2dSrc,
        const long      next2dTrgt,
        const long      current3dSrc,
        const long      current3dTrgt,
        const long      next3dSrc,
        const long      next3dTrgt,
        const long      numProductNodesPerLayer,
        const std::vector< Eigen::Quaterniond>& precRotations) { 

    const double dist2Other2d = getLocalThickness2d(next2dTrgt, features2d, dimFeatures2d);
    const double dist2Other3d = getLocalThickness3d(next3dTrgt, features3d, dimFeatures3d);
    double distDif =  robustLoss(dist2Other2d, dist2Other3d);
    

    double arapL;
    if (PRECOMPUTE_ROTATIONS) {
        arapL =  getArapWPrecompute(numProductNodesPerLayer, currentIdxInLayer,
                nextIdxInLayer, current2dSrc,next2dSrc, precRotations);
    }
    else { 
        std::cout << "PRECOMPUTE_ROTATIONS=False not supported" << std::endl;
    }
    
    return robustLossArap(arapL) + distDif;

}



#endif
