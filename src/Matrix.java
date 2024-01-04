/**
 * {@link Matrix} is a class that can be used for constructing and manipulating real, dense matrices
 * based on the JAMA (Java Matrix) package dated Aug 5, 1998.
 * The original JAMA tools for Singular Value and Eigenvalue Decompositions are included in this same class
 * The class was adapted so it could also be used for sparse matrices, allowing missing entries
 *
 * @author      Antonio Mucherino
 * @author      Simon Hengeveld
 * @since       5 August 2021
 * @version     12 December 2021
 * @see         Maths
 * @package     javaDGP
 */

import java.io.*;
import java.lang.Number;
import java.util.*;
import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

public class Matrix implements Cloneable
{
   /**
    * Row dimension
    */
   private int m;

   /**
    * Column dimension
    */
   private int n;

   /**
    * A two-dimensional array of dimensions m by n containing the values of the {@link Matrix}
    */
   private double [][] A;

   /**
    * A two-dimensional array in which an entry indicates whether the corresponding {@link Matrix#A} entry is set
    */
   private boolean [][] e;

   /**
    * Temporary memory space required for computing the SVD
    */
   private double [][][] SVD;

   /**
    * Temporary memory space required for computing the EVD
    */
   private double [][][] EVD;

   /**
    * Copies an existing {@link Matrix} object in a new one (copy constructor)
    *  @param other         the {@link Matrix} object to copy
    *  @exception IllegalArgumentException {@link Matrix} object is null
    */
   public Matrix(Matrix other)
   {
      if (other == null) throw new IllegalArgumentException("Matrix object is null");
      this.m = other.m;
      this.n = other.n;
      this.e = new boolean [this.m][this.n];
      for (int i = 0; i < this.m; i++) System.arraycopy(other.e[i], 0, this.e[i], 0, this.n);
      this.A = new double [this.m][this.n];
      for (int i = 0; i < this.m; i++) System.arraycopy(other.A[i], 0, this.A[i], 0, this.n);
      this.SVD = null;
      this.EVD = null;
   }

   /**
    * Constructs an empty m-by-n {@link Matrix}
    * @param m          the number of rows
    * @param n          the number of columns
    * @exception IllegalArgumentException Number of rows must be positive
    * @exception IllegalArgumentException Number of columns must be positive
    */
   public Matrix(int m,int n)
   {
      if (m <= 0) throw new IllegalArgumentException("Number of rows must be positive");
      this.m = m;
      if (n <= 0) throw new IllegalArgumentException("Number of columns must be positive");
      this.n = n;
      this.e = new boolean [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.e[i][j] = false;
      this.A = new double [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.A[i][j] = 0.0;
      this.SVD = null;
      this.EVD = null;
   }

   /**
    * Constructs a constant m-by-n {@link Matrix} with constant value
    * @param m          the number of rows
    * @param n          the number of columns
    * @param constant   the constant scalar to fill the {@link Matrix} with
    * @exception IllegalArgumentException Number of rows must be positive
    * @exception IllegalArgumentException Number of columns must be positive
    */
   public Matrix(int m,int n,double constant)
   {
      if (m <= 0) throw new IllegalArgumentException("Number of rows must be positive");
      this.m = m;
      if (n <= 0) throw new IllegalArgumentException("Number of columns must be positive");
      this.n = n;
      this.e = new boolean [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.e[i][j] = true;
      this.A = new double [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.A[i][j] = constant;
      this.SVD = null;
      this.EVD = null;
   }

   /**
    * Constructs a random {@link Matrix} with elements in [0,1]
    * @param m          number of rows
    * @param n          number of columns
    * @param density    (in (0,1]; 1 corresponds to a dense matrix)
    * @exception IllegalArgumentException Density of the {@link Matrix} is not in interval (0,1]
    * @exception IllegalArgumentException Number of rows must be positive
    * @exception IllegalArgumentException Number of columns must be positive
    */
   public Matrix(double density,int m,int n)
   {
      this(m,n);
      if (density <= 0.0 || density > 1.0) throw new IllegalArgumentException("Density of matrix is not in interval (0,1]");

      int k = 0;
      Random R = new Random();
      int N = (int) (density*this.m*this.n);
      if (density > 0.0 && N == 0)  N = 1;
      if (density > 0.5)  N = this.m*this.n - N;
      while (k < N)
      {
         int i = R.nextInt(this.m);
         int j = R.nextInt(this.n);
         if (!this.e[i][j])
         {
            this.e[i][j] = true;
            k++;
         }
      }
      if (density > 0.5)  for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.e[i][j] = !this.e[i][j];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  if (this.e[i][j])  this.A[i][j] = R.nextDouble();
   }

   /**
    * Constructs a diagonal (sparse, square) {@link Matrix} from a 1-D array of primitive doubles
    * @param a          one-dimensional array of doubles
    * @exception IllegalArgumentException The one-dimensional array is null
    */
   public Matrix(double[] a)
   {
      if (a == null) throw new IllegalArgumentException("The one-dimensional array is null");
      this.m = a.length;
      this.n = this.m;
      this.A = new double[m][n];
      this.e = new boolean[m][n];

      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (i == j)
            {
               this.A[i][i] = a[i];
               this.e[i][i] = true;
            }
            else
            {
               this.A[i][j] = 0.0;
               this.e[i][j] = false;
            }
         }
      }
   }

   /**
    * Constructs a diagonal (sparse) {@link Matrix} from a 1-D array of {@link Number} objects
    * @param numb       one-dimensional array of {@link Number} objects
    * @exception IllegalArgumentException the one-dimensional array is null
    */
   public Matrix(Number[] numb)
   {
      if (numb == null) throw new IllegalArgumentException("The one-dimensional array is null");
      this.m = numb.length;
      this.n = this.m;
      this.A = new double[m][n];
      this.e = new boolean[m][n];
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (i == j)
            {
               this.A[i][i] = numb[i].doubleValue();
               this.e[i][i] = true;
            }
            else
            {
               this.A[i][j] = 0.0;
               this.e[i][j] = false;
            }
         }
      }
   }

   /**
    * Constructs a {@link Matrix} from a 2-D array of doubles
    * @param A          two-dimensional array of doubles
    * @exception IllegalArgumentException The two-dimensional array is null
    * @exception IllegalArgumentException The second dimension of the array seems to be zero
    * @exception IllegalArgumentException All rows must have the same length
    */
   public Matrix(double[][] A)
         {
      // verification of the exceptions
      if (A == null) throw new IllegalArgumentException("The two-dimensional array is null");
      this.m = A.length;
      this.n = A[0].length;
      if (this.n == 0) throw new IllegalArgumentException("The second dimension of the array seems to be zero");
      for (int i = 1; i < m; i++)  if (A[i].length != this.n) throw new IllegalArgumentException("All rows must have the same length");

      // loading the data
      this.e = new boolean [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.e[i][j] = true;
      this.A = new double [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.A[i][j] = A[i][j];
      this.SVD = null;
      this.EVD = null;
   }

   /**
    * Construct a {@link Matrix} from a 2-D array of {@link Number} objects
    * @param N          two-dimensional array of {@link Number} objects
    * @exception IllegalArgumentException The two-dimensional array is null
    * @exception IllegalArgumentException The second dimension of the array seems to be zero
    * @exception IllegalArgumentException All rows must have the same length
    */
   public Matrix(Number[][] N)
   {
      // verification of the exceptions
      if (N == null) throw new IllegalArgumentException("The two-dimensional array is null");
      this.m = N.length;
      this.n = N[0].length;
      if (this.n == 0) throw new IllegalArgumentException("The second dimension of the array seems to be zero");
      for (int i = 1; i < m; i++)  if (N[i].length != this.n) throw new IllegalArgumentException("All rows must have the same length");

      // loading the data
      this.e = new boolean [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.e[i][j] = true;
      this.A = new double [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.A[i][j] = N[i][j].doubleValue();
      this.SVD = null;
      this.EVD = null;
   }

   /**
    * Construct a {@link Matrix} from a one-dimensional packed array of doubles
    * @param vals       one-dimensional array of doubles, packed by columns (ala Fortran)
    * @param m          the number of rows
    * @exception IllegalArgumentException Array of doubles is null
    * @exception IllegalArgumentException Number of rows must be positive
    * @exception IllegalArgumentException Array length must be a multiple of m
    */
   public Matrix(double[] vals,int m)
   {
      if (vals == null) throw new IllegalArgumentException("Array of doubles is null");
      int l = vals.length;
      if (m <= 0) throw new IllegalArgumentException("Number of rows must be positive");
      this.m = m;
      this.n = l/m;
      if (l != this.m*this.n) throw new IllegalArgumentException("Array length must be a multiple of m");
      this.e = new boolean [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.e[i][j] = true;
      this.A = new double [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.A[i][j] = vals[i + j*m];
      this.SVD = null;
      this.EVD = null;
   }

   /**
    * Constructs a {@link Matrix} from a one-dimensional packed array of {@link Number} objects
    * @param vals       one-dimensional array of {@link Number} objects, packed by columns (ala Fortran)
    * @param m          number of rows
    * @exception IllegalArgumentException Array of {@link Number} objects is null
    * @exception IllegalArgumentException Number of rows must be positive
    * @exception IllegalArgumentException Array length must be a multiple of m
    */
   public Matrix(Number[] vals,int m)
   {
      if (vals == null) throw new IllegalArgumentException("Array of doubles is null");
      int l = vals.length;
      if (m <= 0) throw new IllegalArgumentException("Number of rows must be positive");
      this.m = m;
      this.n = (m != 0 ? l/m : 0);
      if (l != this.m*this.n) throw new IllegalArgumentException("Array length must be a multiple of m");
      this.e = new boolean [this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.e[i][j] = true;
      this.A = new double[this.m][this.n];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.A[i][j] = vals[i + j*m].doubleValue();
      this.SVD = null;
      this.EVD = null;
   }

   /**
    * Clones the {@link Matrix} object
    * @return The cloned {@link Matrix} object
    */
   public Object clone()
   {
      return new Matrix(this);
   }

   /**
    * Checks if two matrices are approximately equal (with an epsilon)
    * @param other         the  {@link Matrix} to check against
    * @param eps           the tolerance epsilon
    * @exception IllegalArgumentException {@link Matrix} dimensions must agree
    * @return              true if the Matrices are approximately equal, false otherwise
    */
   public boolean epsEqual(Matrix other, double eps){
      checkMatrixDimensions(other);
      for(int i = 0; i < m; i++)
         for(int j = 0; j < n; j++)
            if(!Maths.epsEqual(this.get(i,j), other.get(i,j), eps))
               return false;
      return true;
   }

   /**
    * Simple 'toString' with basic format
    * @return           a String object with the representation of the {@link Matrix}
    */
   public String toString()
   {
      //Note that I use a StringBuilder instead here: concatenating Strings in Java with += is very slow!
      StringBuilder print = new StringBuilder("[");
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (j != 0)  print.append(" ");
            if (this.e[i][j])
            {
               print.append(" ").append(this.A[i][j]);
            }
            else
            {
               print.append(" ?");
            }
         }
         if (i != this.m - 1)
         {
            print.append("\n");
         }
         else
         {
            print.append(" ]");
         }
      }
      return print.toString();
   }

   /**
    * Saves the {@link Matrix} as csv file
    * @param filename   the desired output filename (without extension!)
    */
   public void toCSV(String filename) throws IOException {
      File outputFile = new File(filename + ".csv");
      outputFile.createNewFile();
      FileWriter myWriter = new FileWriter(outputFile);

      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (j != 0)
               myWriter.write(", ");
            if (this.e[i][j])
               myWriter.write(String.valueOf(this.A[i][j]));
            else
               myWriter.write("?");
         }
         if (i != this.m - 1)
            myWriter.write("\n");
      }
      myWriter.close();
   }

   /**
    * Acessor of the number of rows {@link Matrix#m}
    * @return           the number of rows of this {@link Matrix} object
    */
   public int getRowDimension()
   {
      return this.m;
   }

   /**
    * Acessor of the number of columns {@link Matrix#n}
    * @return           the number of columns of this {@link Matrix} object.
    **/
   public int getColumnDimension()
   {
      return this.n;
   }

   /**
    * Returns the number of elements that are set in the {@link Matrix}
    * @return           the number of elements that are set in the {@link Matrix}
    */
   public int getNumberOfElements()
   {
      int count = 0;
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  if (this.e[i][j])  count++;
      return count;
   }

   /**
    * Checks if the {@link Matrix} at hand is a sparse {@link Matrix}
    * @return           true if the {@link Matrix} object is sparse, false otherwise
    */
   public boolean isSparse()
   {
      return this.getNumberOfElements() < this.m*this.n;
   }

   /**
    * Checks if the {@link Matrix} at hand is a dense {@link Matrix}
    * @return           true if the {@link Matrix} object is dense, false otherwise
    */
   public boolean isDense()
   {
      return this.getNumberOfElements() == this.m*this.n;
   }

   /**
    * Checks if the {@link Matrix} at hand is a square {@link Matrix}
    * @return           true if the {@link Matrix} object is square, false otherwise
    */
   public boolean isSquare()
   {
      return this.m == this.n;
   }

   /**
    * Checks if the {@link Matrix} at hand is a hollow {@link Matrix}
    * @return           true if the {@link Matrix} object is hollow, false otherwise
    */
   public boolean isHollow()
   {
      if (!this.isSquare())  return false;
      for (int i = 0; i < this.n; i++)  if (!this.e[i][i] || this.A[i][i] != 0.0)  return false;
      return true;
   }

   /**
    * Checks if the {@link Matrix} at hand is a symmetric {@link Matrix}
    * @return           true if the {@link Matrix} object is symmetric, false otherwise
    */
   public boolean isSymmetric()
   {
      if (!this.isSquare())  return false;
      for (int i = 0; i < this.m; i++)
      {
         for (int j = i + 1; j < this.n; j++)
         {
            if (this.e[i][j] != this.e[j][i])  return false;
            if (this.e[i][j] && this.e[j][i])  if (this.A[i][j] != this.A[j][i])  return false;
         }
      }
      return true;
   }

   /**
    * Make the {@link Matrix} sparse by removing all elements equal to the value in argument
    * @param value      the value of the element to remove everywhere in the {@link Matrix}
    * @exception IllegalStateException The {@link Matrix} is already sparse
    * @return           false if the operation did not succeed (the {@link Matrix} is still dense)
    */
   public boolean sparsify(double value)
   {
      if (this.isSparse()) throw new IllegalStateException("The matrix is already sparse");

      int count = 0;
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (this.A[i][j] == value)
            {
               this.e[i][j] = false;
               count++;
            }
         }
      }
      if (count == 0)  return false;
      return true;
   }

   /**
    * Make the {@link Matrix} dense by setting the missing elements to the value in argument
    * @param value      the constant scalar value for all new elements
    * @exception IllegalStateException The {@link Matrix} is already dense
    **/
   public void densify(double value)
   {
      if (this.isDense()) throw new IllegalStateException("The matrix is already dense");

      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (!this.e[i][j])
            {
               this.A[i][j] = value;
               this.e[i][j] = true;
            }
         }
      }
   }


   /**
    * Checks if a row index is valid (not out of bounds)
    * @param i          the row index to check
    * @return           true if the index is valid, false otherwise
    */
   public boolean isValidRow(int i)
   {
      return i >= 0 && i < this.m;
   }

   /**
    * Checks if a column index is valid (not out of bounds)
    * @param j          the column index to check
    * @return           true if the index is valid, false otherwise
    */
   public boolean isValidColumn(int j)
   {
      return j >= 0 && j < this.n;
   }

   /**
    * Checks if an index pair is valid (not out of bounds)
    * @param i          the row index to check
    * @param j          the column index to check
    * @return           true if the index pair is valid, false otherwise
    */
   public boolean isValidPair(int i,int j)
   {
      return this.isValidRow(i) && this.isValidColumn(j);
   }

   /** Checks if an element of the {@link Matrix} is set
    * @param i          the row index of the element
    * @param j          the column index of the element
    * @exception IllegalArgumentException Indices out of bounds
    @return             true if the element exists, false otherwise
    **/
   public boolean isSet(int i,int j)
   {
      if (!this.isValidPair(i,j)) throw new IllegalArgumentException("Indices out of bounds");
      return this.e[i][j];
   }

   /**
    * Get a {@link Matrix} element
    * @param i          the row index of the element
    * @param j          the column index of the element
    * @exception IllegalArgumentException Indices out of bounds
    * @exception IllegalArgumentException The element (i,j) is not set
    * @return           the element at indices (i,j) as a primitive double
    **/
   public double get(int i,int j)
   {
      if (!this.isValidPair(i,j)) throw new IllegalArgumentException("Indices out of bounds");
      if (!this.e[i][j]) throw new IllegalArgumentException("The element (i,j) is not set");
      return this.A[i][j];
   }

   /**
    * Get the element array {@link Matrix#A} as (primitive) doubles
    * Note: only works if the {@link Matrix} is dense
    * @exception IllegalArgumentException Only for dense matrices
    * @return           a two-dimensional array of unboxed doubles which represent the {@link Matrix}
    */
   public double[][] getArray(){
      if (!this.isDense()) throw new IllegalStateException("Only for dense matrices.");
      return A;

   }

   /**
    * Get the element array {@link Matrix#A} as (wrapper) {@link Double} objects
    * @return           a two-dimensional array of boxed doubles which represent the {@link Matrix}
    */
   public Double[][] getBoxedArray(){
      Double[][] a = new Double[m][n];
      for(int i = 0; i < m; i++)
         for(int j = 0; j < n; j++)
            if(e[i][j])
               a[i][j] = A[i][j];
      return a;
   }

   /**
    * Get the {@link Matrix} row of a given index (only dense version, unknown elements are set to 0.0)
    * @param i          the row index
    * @exception IllegalArgumentException Row index out of bounds
    * @exception IllegalArgumentException Only for dense matrices
    * @return           an array of primitive doubles representing the indicated row
    */
   public double[] getRow(int i)
   {
      if (!this.isValidRow(i)) throw new IllegalArgumentException("Row index out of bounds");
      if (!this.isDense()) throw new IllegalStateException("Only for dense matrices.");

      double [] row = new double [this.n];
      System.arraycopy(this.A[i], 0, row, 0, this.n);
      return row;
   }

   /**
    * Get the {@link Matrix} row of a given index (only dense version, unknown elements are set to null)
    * @param i          the row index
    * @exception IllegalArgumentException Row index out of bounds
    * @return           an array of (wrapper) {@link Double} objects representing the indicated row
    */
   public Double[] getBoxedRow(int i)
   {
      if (!this.isValidRow(i)) throw new IllegalArgumentException("Row index out of bounds");
      Double [] row = new Double [this.n];
      for (int j = 0; j < this.n; j++)
      {
         if (this.e[i][j])  row[j] = this.A[i][j];
      }
      return row;
   }

   /**
    * Get the {@link Matrix} column of a given index (only dense version, unknown elements are set to 0.0)
    * @param j          the column index
    * @exception IllegalArgumentException Column index out of bounds
    * @exception IllegalArgumentException Only for dense matrices
    * @return           an array of primitive doubles representing the indicated column
    */
   public double[] getColumn(int j)
   {
      if (!this.isValidColumn(j)) throw new IllegalArgumentException("Column index out of bounds");
      if (!this.isDense()) throw new IllegalStateException("Only for dense matrices.");

      double [] col = new double [this.m];
      for (int i = 0; i < this.m; i++)
         col[i] = this.A[i][j];
      return col;
   }

   /**
    * Get the {@link Matrix} column of a given index (only dense version, unknown elements are set to null)
    * @param j          the column index
    * @exception IllegalArgumentException Column index out of bounds
    * @return           an array of (wrapper) {@link Double} objects representing the indicated column
    */
   public Double[] getBoxedColumn(int j)
   {
      if (!this.isValidColumn(j)) throw new IllegalArgumentException("Column index out of bounds");
      Double [] col = new Double [this.m];
      for (int i = 0; i < this.m; i++)
      {
         if (this.e[i][j])  col[i] = this.A[i][j];
      }
      return col;
   }

   /**
    * Set the {@link Matrix} element to a given value
    * @param i          the row index of the element
    * @param j          the column index of the element
    * @param value      the new value for the element
    * @exception IllegalArgumentException Indices out of bounds
    * @exception IllegalArgumentException Only a pre-existing element can be set to a new value
    */
   public void set(int i,int j,double value)
   {
      if (!this.isValidPair(i,j)) throw new IllegalArgumentException("Indices out of bounds");
      if (!this.e[i][j]) throw new IllegalArgumentException("Only a pre-existing element can be set to a new value");
      this.A[i][j] = value;
   }

   /**
    * Sets a {@link Matrix} row (only for dense matrices)
    * @param i          the index of the row
    * @param row        the array containing the row
    * @exception IllegalStateException Only for dense matrices
    * @exception IllegalArgumentException The row index is out of bounds
    * @exception IllegalArgumentException The row array is null
    * @exception IllegalArgumentException The row array length does not match with matrix first dimension.
    */
   public void setRow(int i,double[] row)
   {
      // exceptions
      if (!this.isDense()) throw new IllegalStateException("Only for dense matrices");
      if (!this.isValidRow(i)) throw new IllegalArgumentException("The row index is out of bounds");
      if (row == null) throw new IllegalArgumentException("The row array is null");
      if (row.length != this.getColumnDimension()) throw new IllegalArgumentException("The row array length does not match with matrix row length");

      // copy of the row
      for (int j = 0; j < this.n; j++)  this.A[i][j] = row[j];
   }

   /**
    * Sets a {@link Matrix} column (only for dense matrices)
    * @param j         the index of the column
    * @param col       the array containing the column
    * @exception IllegalStateException Only for dense matrices
    * @exception IllegalArgumentException The column index is out of bounds
    * @exception IllegalArgumentException The col array is null
    * @exception IllegalArgumentException The col array length does not match with matrix second dimension
    */
   public void setColumn(int j,double[] col)
   {
      // exceptions
      if (!this.isDense()) throw new IllegalStateException("Only for dense matrices");
      if (!this.isValidColumn(j)) throw new IllegalArgumentException("The column index is out of bounds");
      if (col == null) throw new IllegalArgumentException("The col array is null");
      if (col.length != this.getRowDimension()) throw new IllegalArgumentException("The col array length does not match with matrix second dimension");

      // copy of the column
      for (int i = 0; i < this.m; i++)  this.A[i][j] = col[i];
   }

   /**
    * Adds a new element to a sparse {@link Matrix} object
    * @param i          the row index of the element
    * @param j          the column index of the element
    * @param value      the value of the element
    * @exception IllegalArgumentException Indices out of bounds
    * @exception IllegalArgumentException The element was already set (use 'set' method to change its value)
    */
   public void add(int i,int j,double value)
   {
      if (!this.isValidPair(i,j)) throw new IllegalArgumentException("Indices out of bounds");
      if (this.e[i][j]) throw new IllegalArgumentException("The element was already set (use 'set' method to change its value)");
      this.A[i][j] = value;
      this.e[i][j] = true;
   }

   /**
    * Removes an element from a {@link Matrix} object
    * @param i          the row index of the element
    * @param j          the column index of the element
    * @exception IllegalArgumentException Indices out of bounds
    * @exception IllegalArgumentException The element is not set
    */
   public void remove(int i,int j)
   {
      if (!this.isValidPair(i,j)) throw new IllegalArgumentException("Indices out of bounds");
      if (!this.e[i][j]) throw new IllegalArgumentException("The element is not set");
      this.e[i][j] = false;
   }

   /**
    * Checks whether the {@link Matrix} at hand has the same dimensions as another {@link Matrix}
    * @param other      the other {@link Matrix}
    * @exception IllegalArgumentException {@link Matrix} dimensions must agree
    */
   private void checkMatrixDimensions(Matrix other)
   {
      if (other.m != this.m || other.n != this.n) throw new IllegalArgumentException("Matrix dimensions must agree");
   }

   /**
    * Performs the unary minus operation on the {@link Matrix}
    * @return           the unary minus of the {@link Matrix}
    */
   public Matrix uminus()
   {
      Matrix X = new Matrix(this.m,this.n);
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  X.e[i][j] =  this.e[i][j];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  X.A[i][j] = -this.A[i][j];
      return X;
   }

   /**
    * Computes the transpose of the {@link Matrix}
    * @return           the transpose {@link Matrix}
    */
   public Matrix transpose()
   {
      Matrix X = new Matrix(this.n,this.m);
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  X.e[j][i] = this.e[i][j];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  X.A[j][i] = this.A[i][j];
      return X;
   }

   /**
    * Multiplies the {@link Matrix} by a scalar (C = s * {@link Matrix#A})
    * @param s          the multiplication scalar
    * @return           (C = s * A)
    */
   public Matrix times(double s)
   {
      Matrix X = new Matrix(this.m,this.n);
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  X.e[i][j] = this.e[i][j];
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  X.A[i][j] = s*this.A[i][j];
      return X;
   }

   /**
    * Multiplies the {@link Matrix} by a scalar in place (A = s * {@link Matrix#A})
    * @param s          the multiplication scalar
    * @return           (A = s * A)
    */
   public Matrix timesEquals(double s)
   {
      Matrix X = new Matrix(this.m,this.n);
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  this.A[i][j] = s*this.A[i][j];
      return X;
   }


   /**
    * Computes the trace of the {@link Matrix}
    * @return           the sum of the all the diagonal elements
    */
   public double trace()
   {
      double t = 0.0;
      for (int i = 0; i < Math.min(this.m,this.n); i++)  if (this.e[i][i])  t = t + this.A[i][i];
      return t;
   }

   /**
    * Computes the sum of all the (defined) values of the {@link Matrix}
    * @return           the sum of the all the defined elements
    */
   public double sumValues()
   {
      double sum = 0.0;
      for (int i = 0; i< this.m; i++)  for (int j = 0; j < this.n; j++)  if (this.e[i][j])  sum = sum + this.A[i][j];
      return sum;
   }

   /**
    * Performs a matrix addition (C = {@link Matrix#A} + B)
    * @param B          the other {@link Matrix}
    * @exception IllegalArgumentException {@link Matrix} dimensions must agree
    * @return           (C = A + B)
    */
   public Matrix plus(Matrix B)
   {
      this.checkMatrixDimensions(B);
      Matrix X = new Matrix(m,n);
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  X.e[i][j] = this.e[i][j] || B.e[i][j];
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (X.e[i][j])
            {
               X.A[i][j] = this.A[i][j] + B.A[i][j];
            }
            else if (this.e[i][j])
            {
               X.A[i][j] = this.A[i][j];
            }
            else
            {
               X.A[i][j] = B.A[i][j];
            }
         }
      }
      return X;
   }

   /**
    * Performs a matrix addition in place ({@link Matrix#A} = {@link Matrix#A} + B)
    * @param B          the other {@link Matrix}
    * @exception IllegalArgumentException {@link Matrix} dimensions must agree
    */
   public void plusEquals(Matrix B)
   {
      this.checkMatrixDimensions(B);
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (this.e[i][j] && B.e[i][j])
            {
               this.A[i][j] = this.A[i][j] + B.A[i][j];
            }
            else if (B.e[i][j])
            {
               this.A[i][j] = B.A[i][j];
               this.e[i][j] = true;
            }
         }
      }
   }

   /**
    * Performs a matrix subtraction (C = {@link Matrix#A} - B)
    * @param B          the other {@link Matrix}
    * @exception IllegalArgumentException {@link Matrix} dimensions must agree
    * @return           (C = A - B)
    */
   public Matrix minus(Matrix B)
   {
      this.checkMatrixDimensions(B);
      Matrix X = new Matrix(m,n);
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  X.e[i][j] = this.e[i][j] || B.e[i][j];
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (X.e[i][j])
            {
               X.A[i][j] = this.A[i][j] - B.A[i][j];
            }
            else if (this.e[i][j])
            {
               X.A[i][j] = this.A[i][j];
            }
            else
            {
               X.A[i][j] = -B.A[i][j];
            }
         }
      }
      return X;
   }

   /**
    * Performs a matrix subtraction ({@link Matrix#A} = {@link Matrix#A} - B)
    * @param B          the other {@link Matrix}
    * @exception IllegalArgumentException {@link Matrix} dimensions must agree
    */
   public void minusEquals(Matrix B)
   {
      this.checkMatrixDimensions(B);
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (this.e[i][j] && B.e[i][j])
            {
               this.A[i][j] = this.A[i][j] - B.A[i][j];
            }
            else if (B.e[i][j])
            {
               this.A[i][j] = -B.A[i][j];
               this.e[i][j] = true;
            }
         }
      }
   }

   /**
    * Performs element-by-element multiplication (C= {@link Matrix#A}.*B)
    * @param B          the other {@link Matrix}
    * @exception IllegalArgumentException {@link Matrix} dimensions must agree
    * @return           (C = A.*B)
    */
   public Matrix arrayTimes(Matrix B)
   {
      this.checkMatrixDimensions(B);
      Matrix X = new Matrix(m,n);
      for (int i = 0; i < this.m; i++)  for (int j = 0; j < this.n; j++)  X.e[i][j] = this.e[i][j] && B.e[i][j];
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (X.e[i][j])
            {
               X.A[i][j] = this.A[i][j]*B.A[i][j];
            }
         }
      }
      return X;
   }

   /**
    * Performs element-by-element multiplication ({@link Matrix#A} = {@link Matrix#A}.*B)
    * @param B          the other {@link Matrix}
    * @exception IllegalArgumentException {@link Matrix} dimensions must agree
    */
   public void arrayTimesEquals(Matrix B)
   {
      this.checkMatrixDimensions(B);
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (this.e[i][j] && B.e[i][j])
            {
               this.A[i][j] = this.A[i][j]*B.A[i][j];
            }
            else
            {
               this.e[i][j] = false;
            }
         }
      }
   }

   /**
    * Performs the (Frobernius) inner product (s = sum({@link Matrix#A}.*B))
    * @param B          the other {@link Matrix}
    * @exception IllegalArgumentException {@link Matrix} dimensions must agree
    * @return           (s = sum(A*B))
    */
   public double innerProduct(Matrix B)
   {
      this.checkMatrixDimensions(B);
      Matrix X = this.arrayTimes(B);
      return X.sumValues();
   }

   /**
    * Performs linear algebraic matrix multiplication (A{@link Matrix#A}*B)
    * @param B          the other {@link Matrix}
    * @exception IllegalArgumentException {@link Matrix} inner dimensions must agree
    * @exception IllegalArgumentException {@link Matrix} multiplication allowed only with dense matrices
    * @return           (A * B)
    */
   public Matrix times(Matrix B)
   {
      if (B.m != this.n) throw new IllegalArgumentException("Matrix inner dimensions must agree");
      if (!B.isDense() || !this.isDense()) throw new IllegalArgumentException("Matrix multiplication allowed only with dense matrices");

      Matrix X = new Matrix(this.m,B.n);
      for (int i = 0; i < X.m; i++)  for (int j = 0; j < X.n; j++)  X.e[i][j] = true;

      double[] Bcolj = new double [this.n];
      for (int j = 0; j < B.n; j++)
      {
         for (int k = 0; k < this.n; k++)  Bcolj[k] = B.A[k][j];
         for (int i = 0; i < this.m; i++)
         {
            double[] Arowi = this.A[i];
            double s = 0.0;
            for (int k = 0; k < this.n; k++)  s = s + Arowi[k]*Bcolj[k];
            X.A[i][j] = s;
         }
      }
      return X;
   }

   /**
    * Computes the determinant of the {@link Matrix} (from JAMA LU decomposition)
    * @exception IllegalStateException The determinant is only for dense matrices
    * @exception IllegalStateException The {@link Matrix} must be square
    * @return           the determinant of the {@link Matrix}
    */
   public double determinant () {
      if (this.isSparse()) throw new IllegalStateException("The determinant is only for dense matrices");
      if (m != n) {
         throw new IllegalStateException("The matrix must be square.");
      }
      // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
      double[][] LU = new double[m][n];
      //Using clone does not work well with Double[][] makes shallow copy)
      for (int i = 0; i < m; i++)
         for (int j = 0; j < n; j++)
            LU[i][j] = A[i][j];

      int[] piv = new int[m];
      for (int i = 0; i < m; i++) {
         piv[i] = i;
      }
      int pivsign = 1;
      double[] LUrowi;
      double[] LUcolj = new double[m];

      // Outer loop.
      for (int j = 0; j < n; j++) {
         // Make a copy of the j-th column to localize references.

         for (int i = 0; i < m; i++) {
            LUcolj[i] = LU[i][j];
         }

         // Apply previous transformations.
         for (int i = 0; i < m; i++) {
            LUrowi = LU[i];

            // Most of the time is spent in the following dot product.
            int kmax = Math.min(i,j);
            double s = 0.0;
            for (int k = 0; k < kmax; k++) {
               s += LUrowi[k]*LUcolj[k];
            }

            LUrowi[j] = LUcolj[i] -= s;
         }

         // Find pivot and exchange if necessary.
         int p = j;
         for (int i = j+1; i < m; i++) {
            if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) {
               p = i;
            }
         }
         if (p != j) {
            for (int k = 0; k < n; k++) {
               double t = LU[p][k]; LU[p][k] = LU[j][k]; LU[j][k] = t;
            }
            int k = piv[p]; piv[p] = piv[j]; piv[j] = k;
            pivsign = -pivsign;
         }

         // Compute multipliers.
         if (j < m & LU[j][j] != 0.0) {
            for (int i = j+1; i < m; i++) {
               LU[i][j] /= LU[j][j];
            }
         }
      }

      //The below part comes from the det() method in the LU decomposition code.
      double d = (double) pivsign;
      for (int j = 0; j < n; j++) {
         d *= LU[j][j];
      }
      return d;
   }

   /**
    * Computes the condition number of the {@link Matrix}
    * @exception IllegalStateException Only for dense matrices
    * @return           the condition number computed as max(S)/min(S), where S is the matrix of singular values from SVD.
    */
   public double cond()
   {
      if (this.isSparse()) throw new IllegalStateException("Only for dense matrices");
      double [] s = this.getSingularValues();
      return s[0]/s[Math.min(m,n)-1];
   }

   /**
    *  Computes the effective numerical rank of the {@link Matrix}
    *  @exception IllegalStateException Only for dense matrices
    *  @return             the number of non-negligible singular values
    */
   public int rank()
   {
      if (this.isSparse()) throw new IllegalStateException("Only for dense matrices");
      double [] s = this.getSingularValues();
      double eps = Math.pow(2.0,-52.0);
      double tol = Math.max(m,n)*s[0]*eps;
      int r = 0;
      for (int i = 0; i < s.length; i++)
      {
         if (s[i] > tol)  r++;
      }
      return r;
   }

   /**
    * Computes the 1-norm of the {@link Matrix}
    * @exception IllegalStateException Only for dense matrices
    * @exception IllegalStateException Only for square matrices
    * @return           the maximum of the absolute column sums
    */
   public double norm1()
   {
      if (this.isSparse()) throw new IllegalStateException("Only for dense matrices");
      if (!this.isSquare()) throw new IllegalStateException("Only for square matrices");
      double f = 0.0;
      for (int j = 0; j < this.n; j++)
      {
         double s = 0.0;
         for (int i = 0; i < this.m; i++)
         {
            s += Math.abs(this.A[i][j]);
         }
         f = Math.max(f,s);
      }
      return f;
   }

   /**
    * Computes the 2-norm of the {@link Matrix}
    * @exception IllegalStateException Only for dense matrices
    * @return           the maximal singular value
    */
   public double norm2()
   {
      if (this.isSparse()) throw new IllegalStateException("Only for dense matrices");
      double [] s = this.getSingularValues();
      return s[0];
   }

   /**
    * Computes the infinity-norm of the {@link Matrix}
    * @exception IllegalStateException Only for dense matrices
    * @exception IllegalStateException Only for square matrices
    * @return           the maximum of the absolute row sums
    */
   public double normInf()
   {
      if (this.isSparse()) throw new IllegalStateException("Only for dense matrices");
      if (!this.isSquare()) throw new IllegalStateException("Only for square matrices");
      double f = 0.0;
      for (int i = 0; i < this.m; i++)
      {
         double s = 0.0;
         for (int j = 0; j < this.n; j++)
         {
            s += Math.abs(this.A[i][j]);
         }
         f = Math.max(f,s);
      }
      return f;
   }

   /**
    * Computes the Frobenius-norm of the {@link Matrix}
    * @exception IllegalStateException Only for dense matrices
    * @return           the square root of the sum of squares of all the elements
    */
   public double normF()
   {
      if (this.isSparse()) throw new IllegalStateException("Only for dense matrices");
      double f = 0.0;
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            f = Maths.hypot(f,this.A[i][j]);
         }
      }
      return f;
   }

   /**
    * Computes the singular value decomposition (SVD) of the {@link Matrix} (from JAMA LU decomposition)
    * @exception IllegalStateException SVD is only for dense matrices
    */
   private void SingularValueDecomposition()
   {
      // only for dense matrices
      if (this.isSparse()) throw new IllegalStateException("SVD is only for dense matrices");

      // data
      Matrix A = new Matrix(this);
      int nu = Math.min(this.m,this.n);
      this.SVD = new double [3][][];  // for internal storage; overwrites any previous data
      this.SVD[0] = new double [1][Math.min(this.m+1,this.n)];
      this.SVD[1] = new double [this.m][nu];
      this.SVD[2] = new double [this.n][this.n];

      /* Apparently the failing cases are only a proper subset of (m<n),
	 so let's not throw error.  Correct fix to come later?
      if (m<n) {
	  throw new IllegalArgumentException("Jama SVD only works for m >= n"); }
      */
      double [] s = this.SVD[0][0];
      double [][] U = this.SVD[1];
      double [][] V = this.SVD[2];
      double [] e = new double [this.n];
      double [] work = new double [this.m];
      boolean wantu = true;
      boolean wantv = true;

      // Reduce A to bidiagonal form, storing the diagonal elements
      // in s and the super-diagonal elements in e.
      int nct = Math.min(this.m-1,this.n);
      int nrt = Math.max(0,Math.min(this.n-2,this.m));
      for (int k = 0; k < Math.max(nct,nrt); k++)
      {
         if (k < nct)
         {
            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            s[k] = 0;
            for (int i = k; i < this.m; i++) {
               s[k] = Maths.hypot(s[k],A.A[i][k]);
            }
            if (s[k] != 0.0) {
               if (A.A[k][k] < 0.0) {
                  s[k] = -s[k];
               }
               for (int i = k; i < this.m; i++) {
                  A.A[i][k] /= s[k];
               }
               A.A[k][k] += 1.0;
            }
            s[k] = -s[k];
         }
         for (int j = k+1; j < this.n; j++)
         {
            if ((k < nct) & (s[k] != 0.0))
            {
               // Apply the transformation.
               double t = 0;
               for (int i = k; i < this.m; i++) {
                  t += A.A[i][k]*A.A[i][j];
               }
               t = -t/A.A[k][k];
               for (int i = k; i < this.m; i++) {
                  A.A[i][j] += t*A.A[i][k];
               }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.
            e[j] = A.A[k][j];
         }
         if (wantu & (k < nct)) {

            // Place the transformation in U for subsequent back
            // multiplication.

            for (int i = k; i < this.m; i++) {
               U[i][k] = A.A[i][k];
            }
         }
         if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (int i = k + 1; i < this.n; i++) {
               e[k] = Maths.hypot(e[k],e[i]);
            }
            if (e[k] != 0.0) {
               if (e[k+1] < 0.0) {
                  e[k] = -e[k];
               }
               for (int i = k + 1; i < this.n; i++) {
                  e[i] /= e[k];
               }
               e[k+1] += 1.0;
            }
            e[k] = -e[k];
            if ((k+1 < m) & (e[k] != 0.0))
            {
               // Apply the transformation.
               for (int i = k + 1; i < this.m; i++) {
                  work[i] = 0.0;
               }
               for (int j = k + 1; j < this.n; j++) {
                  for (int i = k + 1; i < this.m; i++) {
                     work[i] += e[j]*A.A[i][j];
                  }
               }
               for (int j = k + 1; j < this.n; j++) {
                  double t = -e[j]/e[k+1];
                  for (int i = k + 1; i < this.m; i++) {
                     A.A[i][j] += t*work[i];
                  }
               }
            }
            if (wantv)
            {
               // Place the transformation in V for subsequent
               // back multiplication.
               for (int i = k + 1; i < this.n; i++) {
                  V[i][k] = e[i];
               }
            }
         }
      }

      // Set up the final bidiagonal matrix or order p.
      int p = Math.min(this.n,this.m+1);
      if (nct < this.n) {
         s[nct] = A.A[nct][nct];
      }
      if (this.m < p) {
         s[p-1] = 0.0;
      }
      if (nrt+1 < p) {
         e[nrt] = A.A[nrt][p-1];
      }
      e[p-1] = 0.0;

      // If required, generate U.
      if (wantu)
      {
         for (int j = nct; j < nu; j++) {
            for (int i = 0; i < this.m; i++) {
               U[i][j] = 0.0;
            }
            U[j][j] = 1.0;
         }
         for (int k = nct-1; k >= 0; k--) {
            if (s[k] != 0.0) {
               for (int j = k + 1; j < nu; j++) {
                  double t = 0.0;
                  for (int i = k; i < this.m; i++) {
                     t += U[i][k]*U[i][j];
                  }
                  t = -t/U[k][k];
                  for (int i = k; i < this.m; i++) {
                     U[i][j] += t*U[i][k];
                  }
               }
               for (int i = k; i < this.m; i++) {
                  U[i][k] = -U[i][k];
               }
               U[k][k] = 1.0 + U[k][k];
               for (int i = 0; i < k - 1; i++) {
                  U[i][k] = 0.0;
               }
            } else {
               for (int i = 0; i < this.m; i++) {
                  U[i][k] = 0.0;
               }
               U[k][k] = 1.0;
            }
         }
      }

      // If required, generate V.
      if (wantv)
      {
         for (int k = this.n - 1; k >= 0; k--)
         {
            if ((k < nrt) & (e[k] != 0.0)) {
               for (int j = k+1; j < nu; j++) {
                  double t = 0;
                  for (int i = k + 1; i < this.n; i++) {
                     t += V[i][k]*V[i][j];
                  }
                  t = -t/V[k+1][k];
                  for (int i = k+1; i < this.n; i++) {
                     V[i][j] += t*V[i][k];
                  }
               }
            }
            for (int i = 0; i < this.n; i++) {
               V[i][k] = 0.0;
            }
            V[k][k] = 1.0;
         }
      }

      // Main iteration loop for the singular values.
      int pp = p-1;
      int iter = 0;
      double eps = Math.pow(2.0,-52.0);
      double tiny = Math.pow(2.0,-966.0);
      while (p > 0)
      {
         int k,kase;

         // Comments coming from JAMA.
         // Here is where a test for too many iterations would go.

         // This section of the program inspects for
         // negligible elements in the s and e arrays.  On
         // completion the variables kase and k are set as follows.

         // kase = 1     if s(p) and e[k-1] are negligible and k<p
         // kase = 2     if s(k) is negligible and k<p
         // kase = 3     if e[k-1] is negligible, k<p, and
         //              s(k), ..., s(p) are not negligible (qr step).
         // kase = 4     if e(p-1) is negligible (convergence).

         for (k = p - 2; k >= -1; k--)
         {
            if (k == -1) {
               break;
            }
            if (Math.abs(e[k]) <=
                    tiny + eps*(Math.abs(s[k]) + Math.abs(s[k+1]))) {
               e[k] = 0.0;
               break;
            }
         }
         if (k == p-2) {
            kase = 4;
         } else {
            int ks;
            for (ks = p-1; ks >= k; ks--) {
               if (ks == k) {
                  break;
               }
               double t = (ks != p ? Math.abs(e[ks]) : 0.) +
                       (ks != k+1 ? Math.abs(e[ks-1]) : 0.);
               if (Math.abs(s[ks]) <= tiny + eps*t)  {
                  s[ks] = 0.0;
                  break;
               }
            }
            if (ks == k) {
               kase = 3;
            } else if (ks == p-1) {
               kase = 1;
            } else {
               kase = 2;
               k = ks;
            }
         }
         k++;

         // Perform the task indicated by kase.
         switch (kase)
         {
            // Deflate negligible s(p).
            case 1: {
               double f = e[p-2];
               e[p-2] = 0.0;
               for (int j = p-2; j >= k; j--) {
                  double t = Maths.hypot(s[j],f);
                  double cs = s[j]/t;
                  double sn = f/t;
                  s[j] = t;
                  if (j != k) {
                     f = -sn*e[j-1];
                     e[j-1] = cs*e[j-1];
                  }
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][p-1];
                        V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
                        V[i][j] = t;
                     }
                  }
               }
            }
            break;

            // Split at negligible s(k).
            case 2: {
               double f = e[k-1];
               e[k-1] = 0.0;
               for (int j = k; j < p; j++) {
                  double t = Maths.hypot(s[j],f);
                  double cs = s[j]/t;
                  double sn = f/t;
                  s[j] = t;
                  f = -sn*e[j];
                  e[j] = cs*e[j];
                  if (wantu) {
                     for (int i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][k-1];
                        U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
                        U[i][j] = t;
                     }
                  }
               }
            }
            break;

            // Perform one qr step.
            case 3: {
               // Calculate the shift.
               double scale = Math.max(Math.max(Math.max(Math.max(
                       Math.abs(s[p-1]),Math.abs(s[p-2])),Math.abs(e[p-2])),
                       Math.abs(s[k])),Math.abs(e[k]));
               double sp = s[p-1]/scale;
               double spm1 = s[p-2]/scale;
               double epm1 = e[p-2]/scale;
               double sk = s[k]/scale;
               double ek = e[k]/scale;
               double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
               double c = (sp*epm1)*(sp*epm1);
               double shift = 0.0;
               if ((b != 0.0) | (c != 0.0)) {
                  shift = Math.sqrt(b*b + c);
                  if (b < 0.0) {
                     shift = -shift;
                  }
                  shift = c/(b + shift);
               }
               double f = (sk + sp)*(sk - sp) + shift;
               double g = sk*ek;

               // Chase zeros.
               for (int j = k; j < p-1; j++) {
                  double t = Maths.hypot(f,g);
                  double cs = f/t;
                  double sn = g/t;
                  if (j != k) {
                     e[j-1] = t;
                  }
                  f = cs*s[j] + sn*e[j];
                  e[j] = cs*e[j] - sn*s[j];
                  g = sn*s[j+1];
                  s[j+1] = cs*s[j+1];
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][j+1];
                        V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
                        V[i][j] = t;
                     }
                  }
                  t = Maths.hypot(f,g);
                  cs = f/t;
                  sn = g/t;
                  s[j] = t;
                  f = cs*e[j] + sn*s[j+1];
                  s[j+1] = -sn*e[j] + cs*s[j+1];
                  g = sn*e[j+1];
                  e[j+1] = cs*e[j+1];
                  if (wantu && (j < m-1)) {
                     for (int i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][j+1];
                        U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
                        U[i][j] = t;
                     }
                  }
               }
               e[p-2] = f;
               iter = iter + 1;
            }
            break;

            // Convergence.
            case 4: {
               // Make the singular values positive.
               if (s[k] <= 0.0) {
                  s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
                  if (wantv) {
                     for (int i = 0; i <= pp; i++) {
                        V[i][k] = -V[i][k];
                     }
                  }
               }

               // Order the singular values.
               while (k < pp) {
                  if (s[k] >= s[k+1]) {
                     break;
                  }
                  double t = s[k];
                  s[k] = s[k+1];
                  s[k+1] = t;
                  if (wantv && (k < n-1)) {
                     for (int i = 0; i < n; i++) {
                        t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
                     }
                  }
                  if (wantu && (k < m-1)) {
                     for (int i = 0; i < m; i++) {
                        t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
                     }
                  }
                  k++;
               }
               iter = 0;
               p--;
            }
            break;
         }
      }
   }

   /**
    * Returns the left singular vectors (from SVD)
    * @exception IllegalStateException SVD is only for dense matrices
    * @return           U
    */
   public Matrix SVD_getU()
   {
      if (this.SVD == null)  this.SingularValueDecomposition();
      return new Matrix(this.SVD[1]);
   }

   /**
    * Returns the right singular vectors (from SVD)
    * @exception IllegalStateException SVD is only for dense matrices
    * @return           V
    */
   public Matrix SVD_getV()
   {
      if (this.SVD == null)  this.SingularValueDecomposition();
      return new Matrix(this.SVD[2]);
   }

   /**
    * Returns the one-dimensional array of (non-trivial) singular values (from SVD)
    * @exception IllegalStateException SVD is only for dense matrices
    * @return           the diagonal array of the singular matrix S
    */
   public double[] getSingularValues()
   {
      if (this.SVD == null)  this.SingularValueDecomposition();
      int N = Math.min(this.m,this.n);
      double [] s = new double [N];
      for (int i = 0; i < N; i++)  s[i] = this.SVD[0][0][i];
      return s;
   }

   /**
    * Returns the diagonal matrix of singular values (from SVD)
    * @exception IllegalStateException SVD is only for dense matrices
    * @return           the diagonal matrix S(igma)
    */
   public Matrix SVD_getS()
   {
      double [] s = this.getSingularValues();
      Matrix S = new Matrix(m,n, 0.0);
      for(int i = 0; i < s.length; i++)
         S.set(i,i, s[i]);
      return S;
   }

   /**
    * Computes the eigenvalue decomposition (EVD) of the {@link Matrix} (from JAMA EVD)
    * @exception IllegalStateException EVD is only for dense matrices
    * @exception IllegalStateException EVD is only for square matrices
    */
   private void EigenvalueDecomposition()
   {
      if (this.isSparse()) throw new IllegalStateException("EVD is only for dense matrices");
      if (!this.isSquare()) throw new IllegalStateException("EVD is only for square matrices");

      // data
      this.EVD = new double [3][][];  // for internal storage; overwrites any previous data
      this.EVD[0] = new double [3][this.n];
      this.EVD[1] = new double [this.n][this.n];
      this.EVD[2] = new double [this.n][this.n];
      double [] d = this.EVD[0][0];
      double [] e = this.EVD[0][1];  // eigenvalues
      double [][] V = this.EVD[1];  // eigenvectors
      double [][] H = this.EVD[2];  // nonsymmetric Hessenberg form
      double [] ort = this.EVD[0][2];  // working space

      // is the matrix symmetric?
      if (this.isSymmetric())
      {
         for (int i = 0; i < this.n; i++)  for (int j = 0; j < this.n; j++)  V[i][j] = this.A[i][j];
         this.tred2();  // Tridiagonalize.
         this.tql2();  // Diagonalize.
      }
      else
      {
         for (int j = 0; j < this.n; j++)  for (int i = 0; i < this.n; i++)  H[i][j] = this.A[i][j];
         this.orthes();  // Reduce to Hessenberg form.
         this.hqr2();  // Reduce Hessenberg to real Schur form.
      }
   }

   /**
    * Computes the tridiagonal form of the {@link Matrix} (used for EVD)
    */
   private void tred2()
   {
      // data
      double [] d = this.EVD[0][0];
      double [] e = this.EVD[0][1];  // eigenvalues
      double [][] V = this.EVD[1];  // eigenvectors
      double [][] H = this.EVD[2];  // nonsymmetric Hessenberg form
      double [] ort = this.EVD[0][2];  // working space

      // Original comment from JAMA package.
      // This is derived from the Algol procedures tred2 by
      // Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
      // Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
      // Fortran subroutine in EISPACK.

      for (int j = 0; j < n; j++) {
         d[j] = V[n-1][j];
      }

      // Householder reduction to tridiagonal form.
      for (int i = n-1; i > 0; i--) {

         // Scale to avoid under/overflow.
         double scale = 0.0;
         double h = 0.0;
         for (int k = 0; k < i; k++) {
            scale = scale + Math.abs(d[k]);
         }
         if (scale == 0.0) {
            e[i] = d[i-1];
            for (int j = 0; j < i; j++) {
               d[j] = V[i-1][j];
               V[i][j] = 0.0;
               V[j][i] = 0.0;
            }
         }
         else
         {
            // Generate Householder vector.
            for (int k = 0; k < i; k++) {
               d[k] /= scale;
               h += d[k] * d[k];
            }
            double f = d[i-1];
            double g = Math.sqrt(h);
            if (f > 0) {
               g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (int j = 0; j < i; j++) {
               e[j] = 0.0;
            }

            // Apply similarity transformation to remaining columns.
            for (int j = 0; j < i; j++) {
               f = d[j];
               V[j][i] = f;
               g = e[j] + V[j][j] * f;
               for (int k = j+1; k <= i-1; k++) {
                  g += V[k][j] * d[k];
                  e[k] += V[k][j] * f;
               }
               e[j] = g;
            }
            f = 0.0;
            for (int j = 0; j < i; j++) {
               e[j] /= h;
               f += e[j] * d[j];
            }
            double hh = f / (h + h);
            for (int j = 0; j < i; j++) {
               e[j] -= hh * d[j];
            }
            for (int j = 0; j < i; j++) {
               f = d[j];
               g = e[j];
               for (int k = j; k <= i-1; k++) {
                  V[k][j] -= (f * e[k] + g * d[k]);
               }
               d[j] = V[i-1][j];
               V[i][j] = 0.0;
            }
         }
         d[i] = h;
      }

      // Accumulate transformations.
      for (int i = 0; i < n-1; i++) {
         V[n-1][i] = V[i][i];
         V[i][i] = 1.0;
         double h = d[i+1];
         if (h != 0.0) {
            for (int k = 0; k <= i; k++) {
               d[k] = V[k][i+1] / h;
            }
            for (int j = 0; j <= i; j++) {
               double g = 0.0;
               for (int k = 0; k <= i; k++) {
                  g += V[k][i+1] * V[k][j];
               }
               for (int k = 0; k <= i; k++) {
                  V[k][j] -= g * d[k];
               }
            }
         }
         for (int k = 0; k <= i; k++) {
            V[k][i+1] = 0.0;
         }
      }
      for (int j = 0; j < n; j++) {
         d[j] = V[n-1][j];
         V[n-1][j] = 0.0;
      }
      V[n-1][n-1] = 1.0;
      e[0] = 0.0;
   }

   /**
    * Computes the tridiagonal form of the {@link Matrix} using a QL method (used for EVD)
    */
   private void tql2()
   {
      // data
      double [] d = this.EVD[0][0];
      double [] e = this.EVD[0][1];  // eigenvalues
      double [][] V = this.EVD[1];  // eigenvectors
      double [][] H = this.EVD[2];  // nonsymmetric Hessenberg form
      double [] ort = this.EVD[0][2];  // working space

      // Original comment from JAMA package.
      // This is derived from the Algol procedures tql2, by
      // Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
      // Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
      // Fortran subroutine in EISPACK.

      for (int i = 1; i < n; i++) {
         e[i-1] = e[i];
      }
      e[n-1] = 0.0;

      double f = 0.0;
      double tst1 = 0.0;
      double eps = Math.pow(2.0,-52.0);
      for (int l = 0; l < n; l++) {

         // Find small subdiagonal element
         tst1 = Math.max(tst1,Math.abs(d[l]) + Math.abs(e[l]));
         int m = l;
         while (m < n) {
            if (Math.abs(e[m]) <= eps*tst1) {
               break;
            }
            m++;
         }

         // If m == l, d[l] is an eigenvalue,
         // otherwise, iterate.
         if (m > l) {
            int iter = 0;
            do {
               iter = iter + 1;  // (Could check iteration count here.)

               // Compute implicit shift
               double g = d[l];
               double p = (d[l+1] - g) / (2.0 * e[l]);
               double r = Maths.hypot(p,1.0);
               if (p < 0) {
                  r = -r;
               }
               d[l] = e[l] / (p + r);
               d[l+1] = e[l] * (p + r);
               double dl1 = d[l+1];
               double h = g - d[l];
               for (int i = l+2; i < n; i++) {
                  d[i] -= h;
               }
               f = f + h;

               // Implicit QL transformation.
               p = d[m];
               double c = 1.0;
               double c2 = c;
               double c3 = c;
               double el1 = e[l+1];
               double s = 0.0;
               double s2 = 0.0;
               for (int i = m-1; i >= l; i--) {
                  c3 = c2;
                  c2 = c;
                  s2 = s;
                  g = c * e[i];
                  h = c * p;
                  r = Maths.hypot(p,e[i]);
                  e[i+1] = s * r;
                  s = e[i] / r;
                  c = p / r;
                  p = c * d[i] - s * g;
                  d[i+1] = h + s * (c * g + s * d[i]);

                  // Accumulate transformation.
                  for (int k = 0; k < n; k++) {
                     h = V[k][i+1];
                     V[k][i+1] = s * V[k][i] + c * h;
                     V[k][i] = c * V[k][i] - s * h;
                  }
               }
               p = -s * s2 * c3 * el1 * e[l] / dl1;
               e[l] = s * p;
               d[l] = c * p;

               // Check for convergence.
            }
            while (Math.abs(e[l]) > eps*tst1);
         }
         d[l] = d[l] + f;
         e[l] = 0.0;
      }

      // Sort eigenvalues and corresponding vectors.
      for (int i = 0; i < n-1; i++) {
         int k = i;
         double p = d[i];
         for (int j = i+1; j < n; j++) {
            if (d[j] < p) {
               k = j;
               p = d[j];
            }
         }
         if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; j++) {
               p = V[j][i];
               V[j][i] = V[j][k];
               V[j][k] = p;
            }
         }
      }
   }

   /**
    * Performs a non-symmetric reduction to Hessenberg form (used for EVD)
    */
   private void orthes()
   {
      // data
      double [] d = this.EVD[0][0];
      double [] e = this.EVD[0][1];  // eigenvalues
      double [][] V = this.EVD[1];  // eigenvectors
      double [][] H = this.EVD[2];  // nonsymmetric Hessenberg form
      double [] ort = this.EVD[0][2];  // working space

      // Original comment from JAMA package.
      // This is derived from the Algol procedures orthes and ortran,
      // by Martin and Wilkinson, Handbook for Auto. Comp.,
      // Vol.ii-Linear Algebra, and the corresponding
      // Fortran subroutines in EISPACK.

      int low = 0;
      int high = n-1;

      for (int m = low+1; m <= high-1; m++) {

         // Scale column.
         double scale = 0.0;
         for (int i = m; i <= high; i++) {
            scale = scale + Math.abs(H[i][m-1]);
         }
         if (scale != 0.0) {

            // Compute Householder transformation.
            double h = 0.0;
            for (int i = high; i >= m; i--) {
               ort[i] = H[i][m-1]/scale;
               h += ort[i] * ort[i];
            }
            double g = Math.sqrt(h);
            if (ort[m] > 0) {
               g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;

            // Apply Householder similarity transformation
            // H = (I-u*u'/h)*H*(I-u*u')/h)
            for (int j = m; j < n; j++) {
               double f = 0.0;
               for (int i = high; i >= m; i--) {
                  f += ort[i]*H[i][j];
               }
               f = f/h;
               for (int i = m; i <= high; i++) {
                  H[i][j] -= f*ort[i];
               }
            }

            for (int i = 0; i <= high; i++) {
               double f = 0.0;
               for (int j = high; j >= m; j--) {
                  f += ort[j]*H[i][j];
               }
               f = f/h;
               for (int j = m; j <= high; j++) {
                  H[i][j] -= f*ort[j];
               }
            }
            ort[m] = scale*ort[m];
            H[m][m-1] = scale*g;
         }
      }

      // Accumulate transformations (Algol's ortran).
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            V[i][j] = (i == j ? 1.0 : 0.0);
         }
      }

      for (int m = high-1; m >= low+1; m--) {
         if (H[m][m-1] != 0.0) {
            for (int i = m+1; i <= high; i++) {
               ort[i] = H[i][m-1];
            }
            for (int j = m; j <= high; j++) {
               double g = 0.0;
               for (int i = m; i <= high; i++) {
                  g += ort[i] * V[i][j];
               }
               // Double division avoids possible underflow
               g = (g / ort[m]) / H[m][m-1];
               for (int i = m; i <= high; i++) {
                  V[i][j] += g * ort[i];
               }
            }
         }
      }
   }

   /**
    * Performs a non-symmetric from Hessenberg to real Schur form (used for EVD)
    */
   private void hqr2()
   {
      // data
      double [] d = this.EVD[0][0];
      double [] e = this.EVD[0][1];  // eigenvalues
      double [][] V = this.EVD[1];  // eigenvectors
      double [][] H = this.EVD[2];  // nonsymmetric Hessenberg form
      double [] ort = this.EVD[0][2];  // working space

      // Original comment from JAMA package.
      // This is derived from the Algol procedure hqr2,
      // by Martin and Wilkinson, Handbook for Auto. Comp.,
      // Vol.ii-Linear Algebra, and the corresponding
      // Fortran subroutine in EISPACK.

      // Initialize
      int nn = this.n;
      int n = nn-1;
      int low = 0;
      int high = nn-1;
      double eps = Math.pow(2.0,-52.0);
      double exshift = 0.0;
      double p=0,q=0,r=0,s=0,z=0,t,w,x,y;

      // Store roots isolated by balanc and compute matrix norm
      double norm = 0.0;
      for (int i = 0; i < nn; i++) {
         if (i < low | i > high) {
            d[i] = H[i][i];
            e[i] = 0.0;
         }
         for (int j = Math.max(i-1,0); j < nn; j++) {
            norm = norm + Math.abs(H[i][j]);
         }
      }

      // Outer loop over eigenvalue index
      int iter = 0;
      while (n >= low) {

         // Look for single small sub-diagonal element
         int l = n;
         while (l > low) {
            s = Math.abs(H[l-1][l-1]) + Math.abs(H[l][l]);
            if (s == 0.0) {
               s = norm;
            }
            if (Math.abs(H[l][l-1]) < eps * s) {
               break;
            }
            l--;
         }

         // Check for convergence
         // One root found
         if (l == n) {
            H[n][n] = H[n][n] + exshift;
            d[n] = H[n][n];
            e[n] = 0.0;
            n--;
            iter = 0;

            // Two roots found
         } else if (l == n-1) {
            w = H[n][n-1] * H[n-1][n];
            p = (H[n-1][n-1] - H[n][n]) / 2.0;
            q = p * p + w;
            z = Math.sqrt(Math.abs(q));
            H[n][n] = H[n][n] + exshift;
            H[n-1][n-1] = H[n-1][n-1] + exshift;
            x = H[n][n];

            // Real pair
            if (q >= 0) {
               if (p >= 0) {
                  z = p + z;
               } else {
                  z = p - z;
               }
               d[n-1] = x + z;
               d[n] = d[n-1];
               if (z != 0.0) {
                  d[n] = x - w / z;
               }
               e[n-1] = 0.0;
               e[n] = 0.0;
               x = H[n][n-1];
               s = Math.abs(x) + Math.abs(z);
               p = x / s;
               q = z / s;
               r = Math.sqrt(p * p+q * q);
               p = p / r;
               q = q / r;

               // Row modification
               for (int j = n-1; j < nn; j++) {
                  z = H[n-1][j];
                  H[n-1][j] = q * z + p * H[n][j];
                  H[n][j] = q * H[n][j] - p * z;
               }

               // Column modification
               for (int i = 0; i <= n; i++) {
                  z = H[i][n-1];
                  H[i][n-1] = q * z + p * H[i][n];
                  H[i][n] = q * H[i][n] - p * z;
               }

               // Accumulate transformations
               for (int i = low; i <= high; i++) {
                  z = V[i][n-1];
                  V[i][n-1] = q * z + p * V[i][n];
                  V[i][n] = q * V[i][n] - p * z;
               }

               // Complex pair
            } else {
               d[n-1] = x + p;
               d[n] = x + p;
               e[n-1] = z;
               e[n] = -z;
            }
            n = n - 2;
            iter = 0;

            // No convergence yet

         } else {

            // Form shift
            x = H[n][n];
            y = 0.0;
            w = 0.0;
            if (l < n) {
               y = H[n-1][n-1];
               w = H[n][n-1] * H[n-1][n];
            }

            // Wilkinson's original ad hoc shift
            if (iter == 10) {
               exshift += x;
               for (int i = low; i <= n; i++) {
                  H[i][i] -= x;
               }
               s = Math.abs(H[n][n-1]) + Math.abs(H[n-1][n-2]);
               x = y = 0.75 * s;
               w = -0.4375 * s * s;
            }

            // MATLAB's new ad hoc shift
            if (iter == 30) {
               s = (y - x) / 2.0;
               s = s * s + w;
               if (s > 0) {
                  s = Math.sqrt(s);
                  if (y < x) {
                     s = -s;
                  }
                  s = x - w / ((y - x) / 2.0 + s);
                  for (int i = low; i <= n; i++) {
                     H[i][i] -= s;
                  }
                  exshift += s;
                  x = y = w = 0.964;
               }
            }

            iter = iter + 1;  // (Could check iteration count here.)

            // Look for two consecutive small sub-diagonal elements
            int m = n-2;
            while (m >= l) {
               z = H[m][m];
               r = x - z;
               s = y - z;
               p = (r * s - w) / H[m+1][m] + H[m][m+1];
               q = H[m+1][m+1] - z - r - s;
               r = H[m+2][m+1];
               s = Math.abs(p) + Math.abs(q) + Math.abs(r);
               p = p / s;
               q = q / s;
               r = r / s;
               if (m == l) {
                  break;
               }
               if (Math.abs(H[m][m-1]) * (Math.abs(q) + Math.abs(r)) <
                       eps * (Math.abs(p) * (Math.abs(H[m-1][m-1]) + Math.abs(z) +
                               Math.abs(H[m+1][m+1])))) {
                  break;
               }
               m--;
            }

            for (int i = m+2; i <= n; i++) {
               H[i][i-2] = 0.0;
               if (i > m+2) {
                  H[i][i-3] = 0.0;
               }
            }

            // Double QR step involving rows l:n and columns m:n
            for (int k = m; k <= n-1; k++) {
               boolean notlast = (k != n-1);
               if (k != m) {
                  p = H[k][k-1];
                  q = H[k+1][k-1];
                  r = (notlast ? H[k+2][k-1] : 0.0);
                  x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                  if (x == 0.0) {
                     continue;
                  }
                  p = p / x;
                  q = q / x;
                  r = r / x;
               }

               s = Math.sqrt(p * p + q * q + r * r);
               if (p < 0) {
                  s = -s;
               }
               if (s != 0) {
                  if (k != m) {
                     H[k][k-1] = -s * x;
                  } else if (l != m) {
                     H[k][k-1] = -H[k][k-1];
                  }
                  p = p + s;
                  x = p / s;
                  y = q / s;
                  z = r / s;
                  q = q / p;
                  r = r / p;

                  // Row modification
                  for (int j = k; j < nn; j++) {
                     p = H[k][j] + q * H[k+1][j];
                     if (notlast) {
                        p = p + r * H[k+2][j];
                        H[k+2][j] = H[k+2][j] - p * z;
                     }
                     H[k][j] = H[k][j] - p * x;
                     H[k+1][j] = H[k+1][j] - p * y;
                  }

                  // Column modification
                  for (int i = 0; i <= Math.min(n,k+3); i++) {
                     p = x * H[i][k] + y * H[i][k+1];
                     if (notlast) {
                        p = p + z * H[i][k+2];
                        H[i][k+2] = H[i][k+2] - p * r;
                     }
                     H[i][k] = H[i][k] - p;
                     H[i][k+1] = H[i][k+1] - p * q;
                  }

                  // Accumulate transformations
                  for (int i = low; i <= high; i++) {
                     p = x * V[i][k] + y * V[i][k+1];
                     if (notlast) {
                        p = p + z * V[i][k+2];
                        V[i][k+2] = V[i][k+2] - p * r;
                     }
                     V[i][k] = V[i][k] - p;
                     V[i][k+1] = V[i][k+1] - p * q;
                  }
               }  // (s != 0)
            }  // k loop
         }  // check convergence
      }  // while (n >= low)

      // Backsubstitute to find vectors of upper triangular form
      if (norm == 0.0) {
         return;
      }

      for (n = nn-1; n >= 0; n--) {
         p = d[n];
         q = e[n];

         // Real vector
         if (q == 0) {
            int l = n;
            H[n][n] = 1.0;
            for (int i = n-1; i >= 0; i--) {
               w = H[i][i] - p;
               r = 0.0;
               for (int j = l; j <= n; j++) {
                  r = r + H[i][j] * H[j][n];
               }
               if (e[i] < 0.0) {
                  z = w;
                  s = r;
               } else {
                  l = i;
                  if (e[i] == 0.0) {
                     if (w != 0.0) {
                        H[i][n] = -r / w;
                     } else {
                        H[i][n] = -r / (eps * norm);
                     }

                     // Solve real equations

                  } else {
                     x = H[i][i+1];
                     y = H[i+1][i];
                     q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                     t = (x * s - z * r) / q;
                     H[i][n] = t;
                     if (Math.abs(x) > Math.abs(z)) {
                        H[i+1][n] = (-r - w * t) / x;
                     } else {
                        H[i+1][n] = (-s - y * t) / z;
                     }
                  }

                  // Overflow control
                  t = Math.abs(H[i][n]);
                  if ((eps * t) * t > 1) {
                     for (int j = i; j <= n; j++) {
                        H[j][n] = H[j][n] / t;
                     }
                  }
               }
            }

            // Complex vector

         } else if (q < 0) {
            int l = n-1;

            // Last vector component imaginary so matrix is triangular

            if (Math.abs(H[n][n-1]) > Math.abs(H[n-1][n])) {
               H[n-1][n-1] = q / H[n][n-1];
               H[n-1][n] = -(H[n][n] - p) / H[n][n-1];
            } else {
               Pair<Double,Double> P = Maths.cdiv(0.0,-H[n-1][n],H[n-1][n-1]-p,q);
               H[n-1][n-1] = (double) P._1();
               H[n-1][n] = (double) P._2();
            }
            H[n][n-1] = 0.0;
            H[n][n] = 1.0;
            for (int i = n-2; i >= 0; i--) {
               double ra,sa,vr,vi;
               ra = 0.0;
               sa = 0.0;
               for (int j = l; j <= n; j++) {
                  ra = ra + H[i][j] * H[j][n-1];
                  sa = sa + H[i][j] * H[j][n];
               }
               w = H[i][i] - p;

               if (e[i] < 0.0) {
                  z = w;
                  r = ra;
                  s = sa;
               } else {
                  l = i;
                  if (e[i] == 0) {
                     Pair<Double,Double> P = Maths.cdiv(-ra,-sa,w,q);
                     H[i][n-1] = (double) P._1();
                     H[i][n] = (double) P._2();
                  } else {

                     // Solve complex equations
                     x = H[i][i+1];
                     y = H[i+1][i];
                     vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                     vi = (d[i] - p) * 2.0 * q;
                     if (vr == 0.0 & vi == 0.0) {
                        vr = eps * norm * (Math.abs(w) + Math.abs(q) +
                                Math.abs(x) + Math.abs(y) + Math.abs(z));
                     }
                     Pair<Double,Double> P = Maths.cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
                     H[i][n-1] = P._1();
                     H[i][n] = P._2();
                     if (Math.abs(x) > (Math.abs(z) + Math.abs(q))) {
                        H[i+1][n-1] = (-ra - w * H[i][n-1] + q * H[i][n]) / x;
                        H[i+1][n] = (-sa - w * H[i][n] - q * H[i][n-1]) / x;
                     } else {
                        P = Maths.cdiv(-r-y*H[i][n-1],-s-y*H[i][n],z,q);
                        H[i+1][n-1] = P._1();
                        H[i+1][n] = P._2();
                     }
                  }

                  // Overflow control
                  t = Math.max(Math.abs(H[i][n-1]),Math.abs(H[i][n]));
                  if ((eps * t) * t > 1) {
                     for (int j = i; j <= n; j++) {
                        H[j][n-1] = H[j][n-1] / t;
                        H[j][n] = H[j][n] / t;
                     }
                  }
               }
            }
         }
      }

      // Vectors of isolated roots
      for (int i = 0; i < nn; i++) {
         if (i < low | i > high) {
            for (int j = i; j < nn; j++) {
               V[i][j] = H[i][j];
            }
         }
      }

      // Back transformation to get eigenvectors of original matrix
      for (int j = nn-1; j >= low; j--) {
         for (int i = low; i <= high; i++) {
            z = 0.0;
            for (int k = low; k <= Math.min(j,high); k++) {
               z = z + V[i][k] * H[k][j];
            }
            V[i][j] = z;
         }
      }
   }

   /**
    * Returns the eigenvector matrix (from EVD)
    * @exception IllegalStateException EVD is only for dense matrices
    * @exception IllegalStateException EVD is only for square matrices
    * @return           V
    */
   public Matrix EVD_getV()
   {
      if (this.EVD == null)  this.EigenvalueDecomposition();
      return new Matrix(this.EVD[1]);
   }

   /**
    * Returns the real parts of the eigenvalues (from EVD)
    * @exception IllegalStateException EVD is only for dense matrices
    * @exception IllegalStateException EVD is only for square matrices.
    * @return           real(diag(D))
    */
   public double[] getRealEigenvalues()
   {
      if (this.EVD == null)  this.EigenvalueDecomposition();
      double [] d = new double [this.n];
      for (int i = 0; i < this.n; i++)  d[i] = this.EVD[0][0][i];
      return d;
   }

   /**
    * Returns the imaginary parts of the eigenvalues (from EVD)
    * @exception IllegalStateException EVD is only for dense matrices
    * @exception IllegalStateException EVD is only for square matrices.
    * @return           imag(diag(D))
    */
   public double[] getImagEigenvalues()
   {
      if (this.EVD == null)  this.EigenvalueDecomposition();
      double [] e = new double [this.n];
      for (int i = 0; i < this.n; i++)  e[i] = this.EVD[0][1][i];
      return e;
   }

   /**
    * Returns the block diagonal eigenvalue matrix (from EVD)
    * @exception IllegalStateException EVD is only for dense matrices
    * @exception IllegalStateException EVD is only for square matrices
    * @return            D
    */
   public Matrix EVD_getD()
   {
      if (this.EVD == null)  this.EigenvalueDecomposition();
      double [] d = this.getRealEigenvalues();
      double [] e = this.getImagEigenvalues();
      Matrix X = Matrix.Zero(this.n,this.n);
      for (int i = 0; i < this.n; i++)
      {
         X.set(i,i,d[i]);
         if (e[i] > 0)  X.A[i][i+1] = e[i];  else if (e[i] < 0)  X.A[i][i-1] = e[i];
      }
      return X;
   }

   /**
    * Prints the {@link Matrix} to stdout, lining the elements up in columns with a Fortran-like 'Fw.d' style format
    * @param w          the column width
    * @param d          the number of digits after the decimal
    */
   public void print (int w, int d) {
      print(new PrintWriter(System.out,true),w,d); }

   /**
    * Prints the {@link Matrix} to an output stream, lines the elements up in columns with a Fortran-like 'Fw.d' style format
    * @param output     the output stream
    * @param w          the column width
    * @param d          the number of digits after the decimal
    */
   public void print (PrintWriter output, int w, int d) {
      DecimalFormat format = new DecimalFormat();
      format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
      format.setMinimumIntegerDigits(1);
      format.setMaximumFractionDigits(d);
      format.setMinimumFractionDigits(d);
      format.setGroupingUsed(false);
      print(output,format,w+2);
   }


   /**
    * Prints the {@link Matrix} to the output stream.  Line the elements up in columns.
    * Uses the format object, and right justify within columns of width
    * characters.
    * Note that is the {@link Matrix} is to be read back in, you probably will want
    * to use a NumberFormat that is set to US Locale.
    * @param output     the output stream
    * @param format     a formatting object to format the {@link Matrix} elements
    * @param width      the column width.
    * @see DecimalFormat#setDecimalFormatSymbols
    */
   public void print (PrintWriter output, NumberFormat format, int width) {
      output.println();  // start on new line.
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            String s = format.format(A[i][j]); // format the number
            int padding = Math.max(1,width-s.length()); // At _least_ 1 space
            for (int k = 0; k < padding; k++)
               output.print(' ');
            output.print(s);
         }
         output.println();
      }
      output.println();   // end with blank line.
   }

   /**
    * Generates the identity matrix (dense)
    * @param m          the number of desired rows
    * @param n          the number of desired columns
    * @exception IllegalArgumentException Number of rows must be positive
    * @exception IllegalArgumentException Number of columns must be positive
    * @return           an m-by-n matrix with ones on the diagonal and zeros elsewhere
    **/
   public static Matrix identity(int m,int n)
   {
      // exceptions
      if (m <= 0) throw new IllegalArgumentException("Number of rows must be positive");
      if (n <= 0) throw new IllegalArgumentException("Number of columns must be positive");

      // constructing the identity matrix (dense version)
      Matrix X = new Matrix(m,n);
      for (int i = 0; i < X.m; i++)
      {
         for (int j = 0; j < X.n; j++)
         {
            X.e[i][j] = true;
         }
      }
      for (int i = 0; i < X.m; i++)
      {
         for (int j = 0; j < X.n; j++)
         {
            X.A[i][j] = (i == j ? 1.0 : 0.0);
         }
      }
      return X;
   }

   /**
    * Generates a zero matrix (dense)
    * @param m          the number of desired rows
    * @param n          the number of desired columns
    * @exception IllegalArgumentException Number of rows must be positive
    * @exception IllegalArgumentException Number of columns must be positive
    * @return           an m-by-n dense matrix with zeros everywhere.
    */
   public static Matrix Zero(int m,int n)
   {
      // exceptions
      if (m <= 0) throw new IllegalArgumentException("Number of rows must be positive");
      if (n <= 0) throw new IllegalArgumentException("Number of columns must be positive");

      // constructing the zero matrix
      Matrix X = new Matrix(m,n);
      for (int i = 0; i < X.m; i++)
      {
         for (int j = 0; j < X.n; j++)
         {
            X.e[i][j] = true;
         }
      }
      for (int i = 0; i < X.m; i++)
      {
         for (int j = 0; j < X.n; j++)
         {
            X.A[i][j] = 0.0;
         }
      }
      return X;
   }

   /**
    * Constructs the EDM from the point-set X (stored as a Matrix object)
    * @param X          the point-set of size d x n
    * @exception IllegalArgumentException X is null
    * @return           the generated EDM
    **/
   public static Matrix EDM(Matrix X)
   {
      int d = X.getRowDimension();
      int n = X.getColumnDimension();

      Matrix gram = X.transpose().times(X);
      double[][] leftExpr = new double[n][n];
      double[][] rightExpr = new double[n][n];

      // Use diagonal of X^T*X
      for (int i = 0; i < n; i++)  for(int j = 0; j < n; j++)  leftExpr[j][i] = rightExpr[i][j] = gram.get(i,i);

      // d_ij = x_i^T * x_i - 2x_i^T*x_j + x_j^T*x_j
      Matrix noisyEDM = new Matrix(leftExpr).minus(gram.times(2).minus(new Matrix(rightExpr)));

      // We make sure that the EDM is symmetric, and we get rid of differences because of rounding
      Matrix M = new Matrix(n,n);
      for (int i = 0; i < n; i++)
      {
         for (int j = i + 1; j < n; j++)
         {
            double val1 = noisyEDM.get(i,j);
            double val2 = noisyEDM.get(i,j);
            double avg = (val1 + val2) / 2;
            M.add(i,j,avg);
            M.set(j,i,avg);  // set entries symmetrically to the average of the two values
         }
      }
      M.densify(0.0);

      return M;
   }
}

