/**
 * {@link DGP} is a class representing a simple undirected distance graph between Cartesian objects
 * We enforce the use of the distance object as weight
 *
 * @author   Antonio Mucherino
 * @author   Simon Hengeveld
 * @since    November 2nd, 2021
 * @version  December, 2023
 * @see      Distance
 * @see      Cartesian
 * @package  ProteinFileReader
 */

import java.util.*;

public class DGP<T extends Cartesian> implements Iterable<T>
{

   /**
    * The list of vertices
    */
   protected ArrayList<T> V;

   /**
    * The edges in a HashMap: pair to edge
    */
   protected HashMap<Pair<T,T>, Distance<T>> E;  // edge set, with potential associated Weight objects

   /**
    * The dimension of the DGP
    */
   protected int K;

   /**
    * Construct the DGP
    * @param K the dimension
    */
   public DGP(int K)
   {
      this.V = new ArrayList<T> ();
      this.E = new HashMap<Pair<T,T>,Distance<T>> ();
      this.K = K;
   }

   /**
    * Returns the dimension of the {@link DGP}
    * @return the dimension
    */
   public int dimension(){
      return this.K;
   }

   /**
    * Counts the number of vertices in the distance graph
    * @return the number of vertices
    */
   public int numberOfVertices()
   {
      return V.size();
   }

   /**
    * Counts the number of edges in the distance graph
    * @return the number of edges
    */
   public int numberOfEdges()
   {
      int nedges = 0;
      int toitself = 0;
      for (Pair<T,T> edge : this.E.keySet())
      {
         T u = edge._1();
         T v = edge._2();
         if (u.equals(v))
            toitself++;
         else
            nedges++;
      }
      //we are an undirected dgp: count every edge only once!
      nedges = nedges/2;

      return toitself + nedges;
   }

   /**
    * Accessor for the list of vertices
    * @return the list of vertices
    */
   public List<T> vertexList()
   {
      ArrayList<T> V = new ArrayList<T>(this.V);
      return V;
   }

   /**
    * Returns the edgeset (containing duplicates, i.e., {u,v} and {v,u} are both included)
    * @return the edgeset
    */
   public Set<Pair<T,T>> edgeSet()
   {
      return this.E.keySet();
   }

   /**
    * Returns a compact edgeset (duplicates are removed, i.e., either {u,v} or {v,u} is included)
    * @return the compacted edgeset
    */
   public Set<Pair<T,T>> edgeSetCompact()
   {
      HashSet<Pair<T,T>> Ec = new HashSet<>();
      for (Pair<T,T> edge : this.E.keySet())
      {
         T u = edge._1();
         T v = edge._2();
         Pair<T,T> other = new Pair<>(v, u);
         if (!Ec.contains(edge) && !Ec.contains(other))  Ec.add(edge);
      }
      return Ec;
   }

   /**
    * Gets the set of distances
    * @return the set of Distances in the graph
    */
   public Set<Distance<T>> distanceSet()
   {
      return new HashSet<>(E.values());
   }


   /**
    * Checks if the graph contains a vertex
    * @param v the vertex
    * @return true if v is contained, false otherwise
    */
   public boolean contains(T v)
   {
      try
      {
         if (v == null) throw new Exception("Vertex given as an argument is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
      return this.V.contains(v);
   }

   /**
    * Checks if the graph contains an edge between u and v
    * @param u the first vertex
    * @param v the second vertex
    * @return true if {u,v} is contained, false otherwise
    */
   public boolean contains(T u,T v)
   {
      try
      {
         if (u == null || v == null) throw new Exception("At least one of the two vertices is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
      Pair<T,T> edge = new Pair<>(u, v);
      return this.E.containsKey(edge);
   }

   public Distance<T> getDistance(Pair<T,T> edge)
   {
      try
      {
         if (edge._1() == null || edge._2() == null) throw new Exception("At least one of the two vertices is null");
         if (!this.contains(edge._1(), edge._2())) throw new Exception("Attempt to access to the weight of an non-existing edge");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
      return this.E.get(edge);

   }

   // weightOf (by two specified vertices)
   public Distance<T> getDistance(T u,T v)
   {
      return getDistance(new Pair<>(u,v));
   }

   /**
    * Get the value of the distance (this throws an exception if it cannot be computed)
    * @param a the first vertex
    * @param b the second vertex
    * @return the reference to the Distance object in the DGP
    */
   public double getExpectedValue(T a, T b){
      return getDistance(a,b).getExpectedValue();
   }


   // addVertex
   // -> the method adds the vertex to the last cluster, 
   //    it creates the cluster0 if DGP is still empty
   public void addVertex(T v)
   {
      try
      {
         if (v == null) throw new Exception("The vertex to be added is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      if (!this.contains(v))
      {
         V.add(v);
      }
   }

   /**
    * Add an edge ({@link Distance}) to the dgp.
    * @param dist the distance to add to the DGP
    */
   public void addEdge(Distance<T> dist){
      T u = dist.getA();
      T v = dist.getB();
      this.E.put(new Pair<>(v,u), dist);
      this.E.put(new Pair<>(u,v), dist);
   }

   // removeEdge (no effect if it doesnt exist)
   public void removeEdge(T u,T v)
   {
      try
      {
         if (u == null || v == null) throw new Exception("At least one of the two vertices is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      Pair<T,T> edge = new Pair<T,T> (u,v);
      this.E.remove(edge);
      edge = new Pair<T,T> (v,u);
      this.E.remove(edge);
   }

   // removeVertex (also removes the corresponding edges, no effect if it doesnt exist)
   public void removeVertex(T v)
   {
      try
      {
         if (v == null) throw new Exception("Vertex given as an argument is null");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      if (this.contains(v))
      {
         // removing the corresponding edges
         ArrayList<Pair<T,T>> toRemove = new ArrayList<>();
         for (Pair<T,T> edge : this.edgeSet())
         {
            T u1 = edge._1();
            T u2 = edge._2();
            if (u1.equals(v) || u2.equals(v))  toRemove.add(edge);
         }
         this.E.keySet().removeAll(toRemove);

         //removing the vertex
         this.V.remove(v);
      }
   }

   /**
    * Checks if the graph is complete
    * @return true if complete, false otherwise
    */
   public boolean isComplete()
   {
      List<T> vertices = this.V;

      for (int i = 0; i < vertices.size(); i++)
         for (int j = i + 1; j < vertices.size(); j++)
            if (!contains(vertices.get(i), vertices.get(j)))  return false;
      return true;
   }

   /**
    * Checks if a vertex order descritizes the DGP search space
    * @param order the vertex order to check
    * @exception IllegalArgumentException The vertex order is attached to a different DGP object
    */
   public boolean isDiscretizable(List<T> order)
   {
      try {
         //we wrap inside a HashSet for better performance!
         if (this.V.size() != order.size() || !new HashSet<>(this.V).containsAll(order)) {
            throw new IllegalArgumentException("The order does not contain all vertices!");
         }
      } catch (Exception e) {
         e.printStackTrace();
         System.exit(1);
      }
      //we need at least K+1 vertices
      if(this.numberOfVertices() < this.K + 1)  return false;

      //first K vertices
      Set<T> firstK = new HashSet<>();
      for(int i = 0; i < this.K; i++)  firstK.add(order.get(i));

      //first K vertices have to be a clique!
      if(!this.subgraph(firstK).isComplete())  return false;

      for(int i = this.K; i < order.size(); i++)
      {
         // we need K reference vertices (u1 < v, u2 < v, ..., uK < v)
         // such that {{u1,v},...,{u2,v},...,{uK,v}} \in E
         int references = 0;
         T v = order.get(i);

         //check for references
         for(int j = i-1; j >= 0; j--){
            T u = order.get(j);
            if(this.contains(u,v)) references++;
            //we are done already
            if(references == this.K) break;
         }

         //not enough reference vertices
         if(references < K)  return false;
      }

      // the instance is discretizable
      return true;
   }

   /**
    * Returns the star of a vertex v in the graph
    * @param v the vertex
    * @return the list of vertices connected to v in the graph
    */
   public List<T> starOf(T v)
   {
      try
      {
         if (v == null) throw new Exception("Specified dgp vertex is null");
         if (!this.contains(v)) throw new Exception("Specified vertex is not contained in the dgp");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }

      ArrayList<T> S = new ArrayList<T> ();
      for (Pair<T,T> edge : this.E.keySet())
      {
         T u1 = edge._1();
         T u2 = edge._2();
         if (v.equals(u1))  if (!S.contains(u2))  S.add(u2);
         if (v.equals(u2))  if (!S.contains(u1))  S.add(u1);
      }
      return S;
   }

   /**
    * The degree of a vertex v in the graph
    * @param v the vertex
    * @return the number of vertices connected to v in the graph
    */
   public int degreeOf(T v)
   {
      return this.starOf(v).size();
   }

   /**
    * Builds a subgraph of the vertices given in a set
    * @param vertexSet the set of vertices
    * @return the subgraph
    */
   public DGP<T> subgraph(Set<T> vertexSet)
   {
      DGP<T> sub = new DGP<T>(this.K);

      // add all vertices with their clusters
      for (T vertex : vertexSet)  sub.addVertex(vertex);

      // add all edges
      for (Pair<T,T> edge : this.edgeSetCompact())
      {
         T u = edge._1();
         T v = edge._2();
         Distance<T> w = this.E.get(edge);
         if (sub.contains(u) && sub.contains(v))  sub.addEdge(w);
      }
      return sub;
   }

   /**
    * Returns an {@link GIterator} attached to this graph
    * @return the iterator
    */
   @Override
   public GIterator iterator()
   {
      return new GIterator(this);
   }

   /**
    * The default iterator for a {@link DGP} instance
    */
   public class GIterator implements Iterator<T>
   {
      /**
       * Reference to the dgp
       */
      protected DGP<T> dgp;

      /**
       * The order
       */
      protected List<T> order;  // the order

      /**
       * The current vertex
       */
      protected int current;

      /**
       * Construct an iterator
       * @param dgp the attached {@link DGP}
       */
      public GIterator(DGP<T> dgp)
      {
         super();
         try
         {
            if (dgp == null) throw new Exception("Impossible to initialize Iterator: the DGP object is null");
            if (V.isEmpty()) throw new Exception("Impossible to initialize Iterator: the DGP is empty");
            this.dgp = dgp;
            this.order = this.dgp.vertexList();
            this.current = -1;
         }
         catch (Exception e)
         {
            e.printStackTrace();
            System.exit(1);
         }
      }

      /**
       * Accessor for the vertex order
       * @return the vertex order
       */
      public List<T> getOrder() {
         return order;
      }
      /**
       * Updates the order to a new one
       * @param order the new order
       * @throws IllegalArgumentException The new order does not contain the same vertices
       *
       */
      public void setOrder(List<T> order){
         //check assumptions and initialize some variables
         try {
            if (order.size() != this.order.size() || !new HashSet<>(order).containsAll(this.order)) {
               throw new IllegalArgumentException("The new order does not contain the same vertices");
            }
         }
         catch (Exception e){
            e.printStackTrace();
            System.exit(1);
         }

         this.order = order;
      }

      /**
       * Check if there is a next vertex in the order
       * @return true if there is a next vertex, false otherwise
       */
      @Override
      public boolean hasNext()
      {
         return this.current + 1 < this.dgp.numberOfVertices();
      }

      /**
       * Accesses the next vertex in the order
       * @return the next vertex in the order
       */
      @Override
      public T next() throws NoSuchElementException
      {
         if (!this.hasNext()) throw new NoSuchElementException();
         this.current++;
         return this.order.get(this.current);
      }

      /**
       * Check if there is a previous vertex in the order
       * @return true if there is a previous vertex, false otherwise
       */
      public boolean hasPrevious()
      {
         return this.current - 1 >= 0;
      }

      /**
       * Accesses the previous vertex in the order
       * @return the previous vertex in the order
       */
      public T previous() throws NoSuchElementException
      {
         if (!this.hasPrevious()) throw new NoSuchElementException();
         this.current--;
         return this.order.get(this.current);
      }

      /**
       * Returns the index of a vertex in the order
       * @param v the vertex to look-up
       * @return the index of the vertex v
       */
      public int indexOf(T v)
      {
         return this.order.indexOf(v);
      }

      /**
       * Moves the iterator to the position of a vertex
       * @param v the vertex to move to
       */
      public void positionAt(T v)
      {
         try
         {
            if (v == null) throw new IllegalArgumentException("The input vertex object is null");
            if (!this.dgp.contains(v)) throw new IllegalArgumentException("The input vertex does not belong to the DGP object of this iterator");
         }
         catch (Exception e)
         {
            e.printStackTrace();
            System.exit(1);
         }

         this.current = this.order.indexOf(v);
      }

      /**
       * Returns a reference to the vertex at an index
       * @param index the index of the vertex
       * @return the vertex at index
       */
      public T get(int index){
         return this.order.get(index);
      }

      /**
       * Returns a sublist of the order
       * @param begin the index first item to be included, inclusive
       * @param end the end index of the sublist, exclusive
       * @return the sublist
       */
      public List<T> subList(int begin, int end){
         return this.order.subList(begin, end);
      }

      /**
       * Resets the iterator
       */
      public void rewind()
      {
         this.current = -1;
      }
   }
}

