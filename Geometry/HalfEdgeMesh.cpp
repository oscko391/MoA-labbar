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
#include "Geometry/HalfEdgeMesh.h"

const unsigned int HalfEdgeMesh::BORDER = (std::numeric_limits<unsigned int>::max)();
const unsigned int HalfEdgeMesh::UNINITIALIZED = (std::numeric_limits<unsigned int>::max)()-1;

HalfEdgeMesh::HalfEdgeMesh()
{
}

HalfEdgeMesh::~HalfEdgeMesh()
{
}

/*! \lab1 Implement the addFace */
/*!
 * \param[in] v1 vertex 1, Vector3<float>
 * \param[in] v2 vertex 2, Vector3<float>
 * \param[in] v3 vertex 3, Vector3<float>
 */
bool HalfEdgeMesh::AddFace(const std::vector<Vector3<float> > &verts){
  // Add your code here
 

  // Add the vertices of the face/triangle
	unsigned int ind0, ind1, ind2;
	AddVertex(verts.at(0), ind0);
	AddVertex(verts.at(1), ind1);
	AddVertex(verts.at(2), ind2);

  // Add all half-edge pairs
	unsigned int indx0, indx1, indx2, indx3, indx4, indx5;
	AddHalfEdgePair(ind0, ind1, indx0, indx1);
	AddHalfEdgePair(ind1, ind2, indx2, indx3);
	AddHalfEdgePair(ind2, ind0, indx4, indx5);
  // Connect inner ring
	e(indx0).next = indx2;
	e(indx0).prev = indx4;
	e(indx2).next = indx4;
	e(indx2).prev = indx0;
	e(indx4).next = indx0;
	e(indx4).prev = indx2;
  // Finally, create the face, don't forget to set the normal (which should be normalized)
	Face f1;
	unsigned int faceInd = mFaces.size(); 

	f1.edge = indx0;
	//calculate normal
	Vector3<float> p1 = v(ind1).pos - v(ind0).pos;
	Vector3<float> p2 = v(ind2).pos - v(ind0).pos;
	
	f1.normal = Cross(p1, p2).Normalize();

	mFaces.push_back(f1);
  // All half-edges share the same left face (previously added)
	e(indx0).face = faceInd;
	e(indx2).face = faceInd;
	e(indx4).face = faceInd;
  // Optionally, track the (outer) boundary half-edges
  // to represent non-closed surfaces
  return true;
}


/*!
 * \param [in] v the vertex to add, Vector3<float>
 * \param [out] indx  the index to the vertex, unsigned int
 * \return a bool indicating whether the HalfEdgeMesh::Vertex was successfully inserted (true) or already existed (false)
 */
bool HalfEdgeMesh::AddVertex(const Vector3<float> & v, unsigned int &indx){
  std::map<Vector3<float>,unsigned int>::iterator it = mUniqueVerts.find(v);
  if (it != mUniqueVerts.end()){
    indx = (*it).second; // get the index of the already existing vertex
    return false;
  }

  mUniqueVerts[v] = indx = GetNumVerts(); // op. [ ] constructs a new entry in map
  Vertex vert;
  vert.pos = v;
  mVerts.push_back(vert); // add it to the vertex list

  return true;
}

/*!
 * Inserts a half edge pair between HalfEdgeMesh::Vertex pointed to by v1 and v2. The first HalfEdgeMesh::HalfEdge (v1->v2) is the inner one, and the second (v2->v1) is the outer.
 * \param [in] v1 vertex 1, Vector3<float>
 * \param [in] v2 vertex 2, Vector3<float>
 * \param [out] indx1  the index to the half-edge from v1 to v2, unsigned int
 * \param [out] indx2  the index to the half-edge from v2 to v1, unsigned int
 * \return a bool indicating whether the half-edge pair was successfully inserted (true) or already existed (false)
 */
bool HalfEdgeMesh::AddHalfEdgePair(unsigned int v1, unsigned int v2, unsigned int &indx1, unsigned int &indx2)
{
   std::map<OrderedPair, unsigned int>::iterator it = mUniqueEdgePairs.find(OrderedPair(v1, v2));
  if (it != mUniqueEdgePairs.end()){
    indx1 = it->second;
    indx2 = e(it->second).pair;
    if(v1 != e(indx1).vert ){
      std::swap(indx1, indx2); // sort correctly
    }
    return false;
  }

  // If not found, calculate both half-edges indices
  indx1 = mEdges.size();
  indx2 = indx1+1;

  // Create edges and set pair index
  HalfEdge edge1, edge2;
  edge1.pair = indx2;
  edge2.pair = indx1;

  // Connect the edges to the verts
  edge1.vert = v1;
  edge2.vert = v2;

  // Connect the verts to the edges
  v(v1).edge = indx1;
  v(v2).edge = indx2;

  // Store the edges in mEdges
  mEdges.push_back(edge1);
  mEdges.push_back(edge2);

  // Store the first edge in the map as an OrderedPair
  OrderedPair op(v1, v2);
  mUniqueEdgePairs[op] = indx1; // op. [ ] constructs a new entry in map, ordering not important
  // sorting done when retrieving

  return true;
}

/*! \lab1 HalfEdgeMesh Implement the MergeAdjacentBoundaryEdge */
/*!
 * Merges the outer UNINITIALIZED/BORDER to an already set inner half-edge.
 * \param [in] indx the index of the INNER half-edge, unsigned int
 */
void HalfEdgeMesh::MergeOuterBoundaryEdge(unsigned int innerEdge)
{
  // Add your code here
  // 1. Merge first loop (around innerEdge->vert)
  // 2. Find leftmost edge, last edge counter clock-wise
  // 3. Test if there's anything to merge
    // 3a. If so merge the gap
    // 3b. And set border flags
  // 4. Merge second loop (around innerEdge->pair->vert)
}

/*! Proceeds to check if the mesh is valid. All indices are inspected and
 * checked to see that they are initialized. The method checks: mEdges, mFaces and mVerts.
 * Also checks to see if all verts have a neighborhood using the findNeighbourFaces method.
 */
void HalfEdgeMesh::Validate()
{
  std::vector<HalfEdge>::iterator iterEdge = mEdges.begin();
  std::vector<HalfEdge>::iterator iterEdgeEnd = mEdges.end();
  while (iterEdge != iterEdgeEnd) {
    if ((*iterEdge).face == UNINITIALIZED ||
        (*iterEdge).next == UNINITIALIZED ||
        (*iterEdge).pair == UNINITIALIZED ||
        (*iterEdge).prev == UNINITIALIZED ||
        (*iterEdge).vert == UNINITIALIZED)
        std::cerr << "HalfEdge " << iterEdge - mEdges.begin() << " not properly initialized" << std::endl;

    iterEdge++;
  }
  std::cerr << "Done with edge check (checked " << GetNumEdges() << " edges)" << std::endl;

  std::vector<Face>::iterator iterTri = mFaces.begin();
  std::vector<Face>::iterator iterTriEnd = mFaces.end();
  while (iterTri != iterTriEnd) {
    if ((*iterTri).edge == UNINITIALIZED)
        std::cerr << "Tri " << iterTri - mFaces.begin() << " not properly initialized" << std::endl;

    iterTri++;
  }
  std::cerr << "Done with face check (checked " << GetNumFaces() << " faces)" << std::endl;

  std::vector<Vertex>::iterator iterVertex = mVerts.begin();
  std::vector<Vertex>::iterator iterVertexEnd = mVerts.end();
  while (iterVertex != iterVertexEnd) {
    if ((*iterVertex).edge == UNINITIALIZED)
        std::cerr << "Vertex " << iterVertex - mVerts.begin() << " not properly initialized" << std::endl;

    iterVertex++;
  }
  std::cerr << "Done with vertex check (checked " << GetNumVerts() << " vertices)" << std::endl;

  std::cerr << "Looping through triangle neighborhood of each vertex... ";
  iterVertex = mVerts.begin();
  iterVertexEnd = mVerts.end();
  int emptyCount = 0;
  std::vector<unsigned int> problemVerts;
  while (iterVertex != iterVertexEnd) {
    std::vector<unsigned int> foundFaces = FindNeighborFaces(iterVertex - mVerts.begin());
    std::vector<unsigned int> foundVerts = FindNeighborVertices(iterVertex - mVerts.begin());
    if (foundFaces.empty() || foundVerts.empty())
      emptyCount++;
    std::set<unsigned int> uniqueFaces(foundFaces.begin(), foundFaces.end());
    std::set<unsigned int> uniqueVerts(foundVerts.begin(), foundVerts.end());
        if ( foundFaces.size() != uniqueFaces.size() ||
         foundVerts.size() != uniqueVerts.size() )
      problemVerts.push_back(iterVertex - mVerts.begin());
    iterVertex++;
  }
  std::cerr << std::endl << "Done: " << emptyCount << " isolated vertices found" << std::endl;
  if(problemVerts.size()){
    std::cerr << std::endl << "Found " << problemVerts.size() << " duplicate faces in vertices: ";
    std::copy(problemVerts.begin(), problemVerts.end(), std::ostream_iterator<unsigned int>(std::cerr, ", "));
    std::cerr << "\n";
  }
  std::cerr << std::endl << "The mesh has genus " << Genus() << ", and consists of " << Shells() << " shells.\n";
}

/*! \lab1 Implement the FindNeighborVertices */
/*! Loops over the neighborhood of a vertex and collects all the vertices sorted counter clockwise.
 * \param [in] vertexIndex  the index to vertex, unsigned int
 * \return a vector containing the indices to all the found vertices.
 */
std::vector<unsigned int> HalfEdgeMesh::FindNeighborVertices(unsigned int vertexIndex) const
{
  // Collected vertices, sorted counter clockwise!
  std::vector<unsigned int> oneRing;
  
  // Add your code here
  Vertex v1 = v(vertexIndex);
  HalfEdge h1 = e(v1.edge);

  //calculate startIndex
  HalfEdge h2 = e(h1.next);
  unsigned int startIndex = h2.vert;

  int counter = 0;
  while (true)
  {
	  h1 = e(h1.next);
	  unsigned int vertInd = h1.vert;

	  if (counter > 0 && vertInd != startIndex)
		  oneRing.push_back(vertInd);
	  else
		  break;

	  h1 = e(h1.next);
	  h1 = e(h1.pair);
	  counter++;
  }
  
  return oneRing;
}

/*! \lab1 Implement the FindNeighborFaces */
/*! Loops over the neighborhood of a vertex and collects all the faces sorted counter clockwise.
 * \param [in] vertexIndex  the index to vertex, unsigned int
 * \return a vector containing the indices to all the found faces.
 */
std::vector<unsigned int> HalfEdgeMesh::FindNeighborFaces(unsigned int vertexIndex) const
{
  // Collected faces, sorted counter clockwise!
  std::vector<unsigned int> foundFaces;

  //find first face
  Vertex v1 = v(vertexIndex);
  HalfEdge h1 = e(v1.edge);
  unsigned int faceInd = h1.face;
  foundFaces.push_back(faceInd);

  unsigned int startIndex= faceInd;
  
  while (true)
  {
	  h1  = e(h1.prev); //gets the halfedges previous pair
	  h1 = e(h1.pair);
	  faceInd = h1.face; //get index to the face of the pair
	  //check if face is already in vector i.e. the ring is complete
	  if (faceInd == startIndex) {
		  break;
	  }
	  foundFaces.push_back(faceInd);
  }

  
  return foundFaces;
}


/*! \lab1 Implement the curvature */
float HalfEdgeMesh::VertexCurvature(unsigned int vertexIndex) const
{
  // Copy code from SimpleMesh or compute more accurate estimate
	std::vector<unsigned int> oneRing = FindNeighborVertices(vertexIndex);
	assert(oneRing.size() != 0);

	unsigned int curr, next;
	const Vector3<float> &vi = mVerts.at(vertexIndex).pos;
	float angleSum = 0;
	float area = 0;
	for (unsigned int i = 0; i<oneRing.size(); i++) {
		// connections
		curr = oneRing.at(i);
		if (i < oneRing.size() - 1)
			next = oneRing.at(i + 1);
		else
			next = oneRing.front();

		// find vertices in 1-ring according to figure 5 in lab text
		// next - beta
		const Vector3<float> &nextPos = mVerts.at(next).pos;
		const Vector3<float> &vj = mVerts.at(curr).pos;

		// compute angle and area
		angleSum += acos((vj - vi)*(nextPos - vi) / ((vj - vi).Length()*(nextPos - vi).Length()));
		area += Cross((vi - vj), (nextPos - vj)).Length()*.5;
	}
	return (2 * M_PI - angleSum) / area;
  
}

float HalfEdgeMesh::FaceCurvature(unsigned int faceIndex) const
{
  // NB Assumes vertex curvature already computed
  unsigned int indx = f(faceIndex).edge;
  const EdgeIterator it = GetEdgeIterator(indx);

  const Vertex& v1 = v(it.GetEdgeVertexIndex());
  const Vertex &v2 = v(it.Next().GetEdgeVertexIndex());
  const Vertex &v3 = v(it.Next().GetEdgeVertexIndex());

  return (v1.curvature + v2.curvature + v3.curvature) / 3.f;
}

Vector3<float> HalfEdgeMesh::FaceNormal(unsigned int faceIndex) const
{
  unsigned int indx = f(faceIndex).edge;
  const EdgeIterator it = GetEdgeIterator(indx);

  const Vector3<float> &p1 = v(it.GetEdgeVertexIndex()).pos;
  const Vector3<float> &p2 = v(it.Next().GetEdgeVertexIndex()).pos;
  const Vector3<float> &p3 = v(it.Next().GetEdgeVertexIndex()).pos;

  const Vector3<float> e1 = p2-p1;
  const Vector3<float> e2 = p3-p1;
  return Cross(e1, e2).Normalize();
}

Vector3<float> HalfEdgeMesh::VertexNormal(unsigned int vertexIndex) const
{

  Vector3<float> n(0,0,0);

  // Add your code here
  std::vector<unsigned int> neighborFaces = FindNeighborFaces(vertexIndex);

 
  for (int i = 0; i < neighborFaces.size(); i++) {
	  n += f(neighborFaces[i]).normal;
  }

  return n.Normalize();
}


void HalfEdgeMesh::Initialize() {
  Validate();
  Update();
}


void HalfEdgeMesh::Update() {
  // Calculate and store all differentials and area

  // First update all face normals and triangle areas
  for(unsigned int i = 0; i < GetNumFaces(); i++){
    f(i).normal = FaceNormal(i);
  }
  // Then update all vertex normals and curvature
  for(unsigned int i = 0; i < GetNumVerts(); i++){
    // Vertex normals are just weighted averages
    mVerts.at(i).normal = VertexNormal(i);
  }

  // Then update vertex curvature
  for(unsigned int i = 0; i < GetNumVerts(); i++){
    mVerts.at(i).curvature = VertexCurvature(i);
    //    std::cerr <<   mVerts.at(i).curvature << "\n";
  }

  // Finally update face curvature
  for(unsigned int i = 0; i < GetNumFaces(); i++){
    f(i).curvature = FaceCurvature(i);
  }

  std::cerr << "Area: " << Area() << ".\n";
  std::cerr << "Volume: " << Volume() << ".\n";

  // Update vertex and face colors
  if (mVisualizationMode == CurvatureVertex) {
    std::vector<Vertex>::iterator iter = mVerts.begin();
    std::vector<Vertex>::iterator iend = mVerts.end();
    float minCurvature = (std::numeric_limits<float>::max)();
    float maxCurvature = -(std::numeric_limits<float>::max)();
    while (iter != iend) {
      if (minCurvature > (*iter).curvature)  minCurvature = (*iter).curvature;
      if (maxCurvature < (*iter).curvature)  maxCurvature = (*iter).curvature;
      iter++;
    }
    std::cerr << "Mapping color based on vertex curvature with range [" << minCurvature << "," << maxCurvature << "]" << std::endl;
    iter = mVerts.begin();
    while (iter != iend) {
      (*iter).color = mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
      iter++;
    }
  }
  else if (mVisualizationMode == CurvatureFace) {
    std::vector<Face>::iterator iter = mFaces.begin();
    std::vector<Face>::iterator iend = mFaces.end();
    float minCurvature = (std::numeric_limits<float>::max)();
    float maxCurvature = -(std::numeric_limits<float>::max)();
    while (iter != iend) {
      if (minCurvature > (*iter).curvature)  minCurvature = (*iter).curvature;
      if (maxCurvature < (*iter).curvature)  maxCurvature = (*iter).curvature;
      iter++;
    }
    std::cerr << "Mapping color based on face curvature with range [" << minCurvature << "," << maxCurvature << "]" << std::endl;
    iter = mFaces.begin();
    while (iter != iend) {
      (*iter).color = mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
      iter++;
    }
  }

}


/*! \lab1 Implement the area */
float HalfEdgeMesh::Area() const
{
  float area = 0;
  // Add code here
  //loop over all faces
  for (int i = 0; i < mFaces.size(); i++)
  {
	  HalfEdge h1 = e(f(i).edge); //get the faces edge
	  Vector3<float> v1 = v(h1.vert).pos;

	  HalfEdge h2 = e(h1.next);
	  Vector3<float> v2 = v(h2.vert).pos;

	  HalfEdge h3 = e(h2.next);
	  Vector3<float> v3 = v(h3.vert).pos;

	  Vector3<float> p1 = v2 - v1;
	  Vector3<float> p2 = v3 - v1;
		
	  Vector3<float> crossProd = Cross(p1, p2);
	  area += crossProd.Length() / 2;
  }
  std::cout << "Area : " << area << std::endl;
  return area;
}

/*! \lab1 Implement the volume */
float HalfEdgeMesh::Volume() const
{
  float volume = 0;
  // Add code here
  for (int i = 0; i < mFaces.size(); i++)
  {
	  HalfEdge h1 = e(f(i).edge); //get the faces edge
	  Vector3<float> v1 = v(h1.vert).pos;
	  
	  HalfEdge h2 = e(h1.next);
	  Vector3<float> v2 = v(h2.vert).pos;

	  HalfEdge h3 = e(h2.next);
	  Vector3<float> v3 = v(h3.vert).pos;

	  Vector3<float> midPoint = (v1 + v2 + v3) / 3;
	  //volume += midPoint * f(h1.face).normal;

	  Vector3<float> p1 = v2 - v1;
	  Vector3<float> p2 = v3 - v1;

	  Vector3<float> crossProd = Cross(p1, p2);
	  float area = crossProd.Length() / 2;

	  volume += midPoint * f(h1.face).normal * area;
  }

  
  return volume/3;
}

/*! \lab1 Calculate the number of shells  */
int HalfEdgeMesh::Shells() const
{
  return 1;
}

/*! \lab1 Implement the genus */
int HalfEdgeMesh::Genus() const
{
  // Add code here
  std::cerr << "Genus calculation not implemented for half-edge mesh!\n";
  return 0;
}

void HalfEdgeMesh::Dilate(float amount)
{
  std::vector<Vertex>::iterator iter = mVerts.begin();
  std::vector<Vertex>::iterator iend = mVerts.end();
  while (iter != iend) {
    (*iter).pos += amount*(*iter).normal;
    iter++;
  }

  Initialize();
  Update();
}

void HalfEdgeMesh::Erode(float amount)
{
  std::vector<Vertex>::iterator iter = mVerts.begin();
  std::vector<Vertex>::iterator iend = mVerts.end();
  while (iter != iend) {
    (*iter).pos -= amount*(*iter).normal;
    iter++;
  }

  Initialize();
  Update();
}

void HalfEdgeMesh::Smooth(float amount)
{
  std::vector<Vertex>::iterator iter = mVerts.begin();
  std::vector<Vertex>::iterator iend = mVerts.end();
  while (iter != iend) {
    (*iter).pos -= amount*(*iter).normal * (*iter).curvature;
    iter++;
  }

  Initialize();
  Update();
}

void HalfEdgeMesh::Render()
{
  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  // Apply transform
  glPushMatrix(); // Push modelview matrix onto stack

  // Convert transform-matrix to format matching GL matrix format
  // Load transform into modelview matrix
  glMultMatrixf(  mTransform.ToGLMatrix().GetArrayPtr() );

  if (mWireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // Draw geometry
  glBegin(GL_TRIANGLES);
  const int numTriangles = GetNumFaces();
  for (int i = 0; i < numTriangles; i++){

    Face & face = f(i);

    HalfEdge* edge = &e(face.edge);

    Vertex & v1 = v(edge->vert);
    edge = &e(edge->next);

    Vertex & v2 = v(edge->vert);
    edge = &e(edge->next);

    Vertex & v3 = v(edge->vert);

    if (mVisualizationMode == CurvatureVertex) {
      glColor3fv(v1.color.GetArrayPtr());
      glNormal3fv(v1.normal.GetArrayPtr());
      glVertex3fv(v1.pos.GetArrayPtr());

      glColor3fv(v2.color.GetArrayPtr());
      glNormal3fv(v2.normal.GetArrayPtr());
      glVertex3fv(v2.pos.GetArrayPtr());

      glColor3fv(v3.color.GetArrayPtr());
      glNormal3fv(v3.normal.GetArrayPtr());
      glVertex3fv(v3.pos.GetArrayPtr());
    }
    else {
      glColor3fv(face.color.GetArrayPtr());
      glNormal3fv(face.normal.GetArrayPtr());

      glVertex3fv(v1.pos.GetArrayPtr());
      glVertex3fv(v2.pos.GetArrayPtr());
      glVertex3fv(v3.pos.GetArrayPtr());
    }

  }
  glEnd();

  // Mesh normals by courtesy of Richard Khoury
  if (mShowNormals)
  {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    const int numTriangles = GetNumFaces();
    for (int i = 0; i < numTriangles; i++){

      Face & face = f(i);

      HalfEdge* edge = &e(face.edge);

      Vertex & v1 = v(edge->vert);
      edge = &e(edge->next);

      Vertex & v2 = v(edge->vert);
      edge = &e(edge->next);

      Vertex & v3 = v(edge->vert);

      Vector3<float> faceStart = (v1.pos + v2.pos + v3.pos) / 3.0f;
      Vector3<float> faceEnd = faceStart + face.normal*0.1;

      glColor3f(1,0,0); // Red for face normal
      glVertex3fv(faceStart.GetArrayPtr());
      glVertex3fv(faceEnd.GetArrayPtr());

      glColor3f(0,1,0); // Vertex normals in Green
      glVertex3fv(v1.pos.GetArrayPtr());
      glVertex3fv((v1.pos + v1.normal*0.1).GetArrayPtr());
      glVertex3fv(v2.pos.GetArrayPtr());
      glVertex3fv((v2.pos + v2.normal*0.1).GetArrayPtr());
      glVertex3fv(v3.pos.GetArrayPtr());
      glVertex3fv((v3.pos + v3.normal*0.1).GetArrayPtr());
    }
    glEnd();
    glEnable(GL_LIGHTING);

  }


  if (mWireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // Restore modelview matrix
  glPopMatrix();

  GLObject::Render();

}

