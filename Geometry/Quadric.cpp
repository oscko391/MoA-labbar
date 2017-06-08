/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "Quadric.h"

Quadric::Quadric(const Matrix4x4<float> & q){
  this->mQuadric = q;
}

Quadric::~Quadric(){}

/*!
 * Evaluation of world coordinates are done through either transformation
 * of the world-coordinates by mWorld2Obj, or transformation of the quadric
 * coefficient matrix by GetTransform() ONCE (see Section 2.2 in lab text).
 */
float Quadric::GetValue(float x, float y, float z) const
{
	Implicit::TransformW2O(x, y, z);
	Vector4<float> position4D = Vector4<float> (x, y, z, 1.0);


	// p^tQp = f(x,y,z)
  return position4D * (mQuadric * position4D);
}

/*!
 * Use the quadric matrix to evaluate the gradient.
 */
Vector3<float> Quadric::GetGradient(float x, float y, float z) const
{
	Implicit::TransformW2O(x, y, z);
	Vector4<float> position4D = Vector4<float>(x, y, z, 1.0);
	Vector4<float> result = 2.0f * (mQuadric*position4D);

  return Vector3<float>(result[0], result[1], result[2]);
}

