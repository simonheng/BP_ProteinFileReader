/**
 * {@link Revorder} is a class extending {@link DGP.GIterator}, representing an iterator for an order that discretizes the search space of the DGP.
 * See: Jérémy Omer; Antonio Mucherino. The Referenced Vertex Ordering Problem: Theory, Applications, and Solution Methods. Open Journal of Mathematical Optimization, Volume 2 (2021), article no. 6
 * @author  Antonio Mucherino
 * @author  Simon Hengeveld
 * @since   2021
 * @version December 2023
 * @see     Cartesian
 * @see     Coordinates
 * @package ProteinFileReader
 */

import java.util.*;

public class Revorder<T extends Cartesian> extends DGP<T>.GIterator
{
   /**
    * The name of the current order ("original","greedy")
    */
   private String name;

   /**
    * The {@link DGP} dimension
    */
   private int K = 3;

   /**
    * The current Kplet
    */
   private int[] kplet;

   /**
    * Constructs a Revorder: if the current order (in dgp.vertexList()) is not discretizable, use the greedy method to attempt to discretize
    * @param dgp the attached {@link DGP} object
    */
   public Revorder(DGP<T> dgp)
   {
      dgp.super(dgp);
      try
      {
         this.kplet = null;
         if (dgp.isDiscretizable(dgp.vertexList()))  // CAREFUL: isDiscretizable does not distinguish between exact and interval distances!
            this.name = "original";
         else if (!this.greedy())
            throw new IllegalStateException("Revorder does not exist for the given DGP instance");
         this.rewind();
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * Constructs a Revorder using a proposed order. If the proposal is not discretizable, use the greedy method to attempt to discretize
    * @param dgp the attached {@link DGP} object
    * @param order the proposed order
    */
   public Revorder(DGP<T> dgp, List<T> order)
   {
      dgp.super(dgp);
      try
      {
         this.K = dgp.dimension();
         this.kplet = null;
         this.setOrder(order);
         if (dgp.isDiscretizable(order))  // CAREFUL: isDiscretizable does not distinguish between exact and interval distances!
            this.name = "original";
         else if (!this.greedy())
            throw new IllegalStateException("Revorder does not exist for the given DGP instance");
         this.rewind();
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
   }

   /**
    * Resets the iterator
    */
   @Override
   public void rewind()
   {
      super.rewind();
      this.kplet = null;
   }

   /**
    * Reorder the vertices
    * @param method the reordering method (only greedy supported for now)
    * @return true if the reordering was succesful, false otherwise
    */
   public boolean reorder(String method)
   {
      try
      {
         if (method == null) throw new Exception("String argument supposed to contain the ordering method is null");
         if (method.equalsIgnoreCase("greedy"))
            return this.greedy();
         else
            throw new IllegalArgumentException("Specified ordering method is unknown");
      }
      catch (Exception e)
      {
         e.printStackTrace();
         System.exit(1);
      }
      return false;
   }

   /**
    * Reorder the vertices using the greedy approach
    * @return true if the reordering was succesful, false otherwise
    */
   private boolean greedy()
   {
      // using current Revorder to go over current order and identify the K-cliques in the dgp
      this.rewind();
      if (!this.hasNextKplet())  return false;

      // trying the first identified clique
      DGP<T> clique = this.nextClique(true);
      System.out.println(clique);

      if (clique == null)  return false;


      do  // looping over the other cliques if the first one doesnt lead to any valid order
      {
         List<T> cliqueVertexSet = clique.vertexList();
         ArrayList<T> newOrder = new ArrayList<> (this.dgp.numberOfVertices());
         newOrder.addAll(cliqueVertexSet);  // the order of the vertices in the clique is not important
         Set<T> remaining = new HashSet<T> (this.order);
         remaining.removeAll(cliqueVertexSet);

         // constructing the rest of the ordering
         while (remaining.size() > 0)
         {
            T max = null;
            int maxRefs = K - 1;
            for (T v : this.order)
            {
               if (!newOrder.contains(v))
               {
                  Set<T> neighbours = new HashSet<>(this.dgp.starOf(v));
                  neighbours.removeAll(remaining);
                  int nRefs = neighbours.size();
                  if (nRefs > maxRefs)
                  {
                     int neRefs = 0;
                     for (T u : neighbours)
                     {
                        Distance<T> d = this.dgp.getDistance(u,v);
                        if (d.hasExpectedValue())  neRefs++;
                     }
                     if (neRefs >= K - 1)
                     {
                        maxRefs = nRefs;
                        max = v;
                     }
                  }
               }
            }
            if (max != null)
            {
               newOrder.add(max);
               remaining.remove(max);
            }
            else break;
         }

         // were we able to construct the ordering?
         if (newOrder.size() == this.dgp.numberOfVertices())
         {
            this.order = newOrder;
            this.name = "greedy";
            return true;
         }
         else
         {
            // trying to detect another clique
            clique = this.nextClique(true);
         }
      }
      while (clique != null);

      // ending
      return false;
   }

   /**
    * Check if there is a nextKplet
    * @return true if there is a nextKplet, false otherwise
    */
   public boolean hasNextKplet()
   {
      if (this.kplet == null)  return true;
      for (int i = 1; i <= this.K; i++)
      {
         if (this.kplet[this.K - i] < this.dgp.numberOfVertices() - i)  return true;
      }
      return false;
   }

   /**
    * Accesses the next Kplet in the order
    * @return the next Kplet
    */
   public DGP<T> nextKplet()
   {
      // is this the first kplet?
      boolean done = false;
      if (this.kplet == null)
      {
         this.kplet = new int [this.K];
         for (int k = 0; k < this.K; k++)  this.kplet[k] = k;
         done = true;
      }
      else
      {
         // defining the next kplet
         int i = 1;
         while (i <= this.K && !done)
         {
            int k = this.K - i;
            if (this.kplet[k] < this.dgp.numberOfVertices() - i)
            {
               this.kplet[k]++;
               for (int h = k + 1; h < this.K; h++)  this.kplet[h] = this.kplet[h-1] + 1;
               done = true;
            }
            i++;
         }
      }

      // does the kplet really exist?
      if (!done)  return null;

      // constructing the subgraph corresponding to the kplet
      Set<T> vertexSet = new HashSet<> ();
      for (int i = 0; i < this.K; i++)  vertexSet.add(this.order.get(this.kplet[i]));
      return this.dgp.subgraph(vertexSet);
   }

   /**
    * Accesses the next clique
    * @param allExact whether we want to require all distances to be exact
    * @return the next clique
    */
   public DGP<T> nextClique(boolean allExact)
   {
      while (this.hasNextKplet())
      {
         DGP<T> sgraph = this.nextKplet();
         if (sgraph != null && sgraph.isComplete())
         {
            // if we don't need to have exact distances, we are done
            if (!allExact)  return sgraph;

            // otherwise, check if all involved distances are exact
            boolean areExact = true;
            for (Pair<T,T> edge : sgraph.edgeSetCompact())
            {
               Distance<T> dist = sgraph.getDistance(edge);
               if (!dist.hasExpectedValue())
               {
                  areExact = false;
                  break;
               }
            }
            if (areExact)  return sgraph;
         }
      }
      return null;
   }

   /**
    * Accesses the next (non-exact) clique
    * @return the next (non-exact) clique
    */
   public DGP<T> nextClique()
   {
      return this.nextClique(false);
   }

}

