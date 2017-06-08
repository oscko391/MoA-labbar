/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "QuadricDecimationMesh.h"

const QuadricDecimationMesh::VisualizationMode QuadricDecimationMesh::QuadricIsoSurfaces = NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize()
{
  // Allocate memory for the quadric array
  unsigned int numVerts = mVerts.size();
  mQuadrics.reserve(numVerts);
  std::streamsize width = std::cerr.precision(); // store stream precision
  for (unsigned int i = 0; i < numVerts; i++) {

    // Compute quadric for vertex i here
    mQuadrics.push_back(createQuadricForVert(i));


    // Calculate initial error, should be numerically close to 0

    Vector3<float> v0 = mVerts[i].pos;
    Vector4<float> v(v0[0],v0[1],v0[2],1);
    Matrix4x4<float> m = mQuadrics.back();

    float error = v*(m*v);
    //std::cerr << std::scientific << std::setprecision(2) << error << " ";
  }
  std::cerr << std::setprecision(width) << std::fixed; // reset stream precision

  // Run the initialize for the parent class to initialize the edge collapses
  DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute, DecimationMesh::EdgeCollapse
 */
void QuadricDecimationMesh::computeCollapse(EdgeCollapse * collapse)
{
	
  // Compute collapse->position and collapse->cost here
  // based on the quadrics at the edge endpoints
	unsigned int vertIndex1 = e(collapse->halfEdge).vert;
	unsigned int vertIndex2 = e(e(collapse->halfEdge).pair).vert;
	
	Matrix4x4<float> Q1 = mQuadrics[vertIndex1];
	Matrix4x4<float> Q2 = mQuadrics[vertIndex2];
	Matrix4x4<float> Q = Q1 + Q2;
	Matrix4x4<float> Qho = Q;


	Qho(3, 0) = Qho(3, 1) = Qho(3, 2) = 0.0f;
	Qho(3, 3) = 1.0f;
	
	//calculate cost and pos if matrix is invertable
	if (!Qho.IsSingular()) {
		Vector4<float> vec(0.0f, 0.0f, 0.0f, 1.0f);
		//compute new position
		vec = Qho.Inverse()*vec;
		collapse->position = Vector3<float>(vec[0],vec[1],vec[2]);

		//compute cost
		float deltaV = vec*(Q*vec);
		collapse->cost = deltaV;
	}
	//calculate cost a pos for non invertable matrix
	else {
		//v1
		Vector3<float> vec1 = v(vertIndex1).pos;
		Vector4<float> v1( vec1[0],vec1[1],vec1[2], 1.0f);
		float v1Cost = v1*(Q*v1);
		//v2
		Vector3<float> vec2 = v(vertIndex2).pos;
		Vector4<float> v2(vec2[0], vec2[1], vec2[2], 1.0f);
		float v2Cost = v2*(Q*v2);
		//(v1+v2)/2

		Vector3<float> midVec = (v(vertIndex1).pos + v(vertIndex2).pos) / 2.0f;
		Vector4<float> midV(midVec[0], midVec[1], midVec[2], 1.0f);
		float v3Cost = midV*(Q*midV);

		//collapse->cost = std::min(std::min(v1Cost, v2Cost), v3Cost);
		collapse->position = vec1;
		collapse->cost = v1Cost;
		if (v2Cost < v1Cost) {
			collapse->position = vec2;
			collapse->cost = v2Cost;
		}
		if (v3Cost < collapse->cost) {
			collapse->position = midVec;
			collapse->cost = v3Cost;
		}

	}
	
	
	
  
}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(unsigned int ind)
{
  DecimationMesh::updateVertexProperties(ind);
  mQuadrics[ind] = createQuadricForVert(ind);
}

/*!
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
Matrix4x4<float> QuadricDecimationMesh::createQuadricForVert(unsigned int indx) const{
  float q[4][4] = {{0,0,0,0},
                   {0,0,0,0},
                   {0,0,0,0},
                   {0,0,0,0}};
  Matrix4x4<float> Q(q);

  // The quadric for a vertex is the sum of all the quadrics for the adjacent faces

  //loop over neighbour faces and add to quadric sum
  std::vector<unsigned int> neighborFaces = FindNeighborFaces(indx);

  for (int i = 0; i < neighborFaces.size(); i++)
  {
	  Q += createQuadricForFace(neighborFaces[i]);
  }

  // Tip: Matrix4x4 has an operator +=
  return Q;
}

/*!
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
Matrix4x4<float> QuadricDecimationMesh::createQuadricForFace(unsigned int indx) const{

  // Calculate the quadric (outer product of plane parameters) for a face
  // here using the formula from Garland and Heckbert
	//a,b and c of the planes equation
	Vector3<float> planeNorm = f(indx).normal;

	//x,y,z for a point on the plane
	Vector3<float> position = v(e(f(indx).edge).vert).pos;
	//d = -(ax+by+cz)
	float d = -(position[0] * planeNorm[0] + position[1] * planeNorm[1] + position[2] * planeNorm[2]);

	float K[4][4] = 
	{
		{planeNorm[0] * planeNorm[0],planeNorm[0] * planeNorm[1], planeNorm[0] * planeNorm[2], planeNorm[0] * d },
		{ planeNorm[0] * planeNorm[1], planeNorm[1] * planeNorm[1], planeNorm[1] * planeNorm[2], planeNorm[1] * d },
		{ planeNorm[0] * planeNorm[2], planeNorm[1] * planeNorm[2], planeNorm[2] * planeNorm[2], planeNorm[2] * d },
		{ planeNorm[0] *d, planeNorm[1] *d, planeNorm[2] * d, d * d }
	};


  return Matrix4x4<float>(K);
}


void QuadricDecimationMesh::Render()
{
  DecimationMesh::Render();

  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  if (mVisualizationMode == QuadricIsoSurfaces)
    {
      // Apply transform
      glPushMatrix(); // Push modelview matrix onto stack

      // Implement the quadric visualization here
      std::cout << "Quadric visualization not implemented" << std::endl;

      // Restore modelview matrix
      glPopMatrix();
    }
}

