
/**
 * {@link Geometry} is a static class containing methods that operate on geometrical objects
 *
 * @author      Antonio Mucherino
 * @author      Simon Hengeveld
 * @since       November 2nd, 2021
 * @version     December, 2023
 * @see         Maths
 * @package     ProteinFileReader
 */


public class Geometry
{
    /**
     * Computes the centroid of a point-set
     * @param points    the point-set
     * @return          the centroid (coordinates represented by an array)
     */
    public static double[] centroid(double[][] points) {
        double[] centroid = new double[points[0].length];

        //We sum the coordinates of all the points for each of the coordinate axes
        for (double[] point : points) {
            for (int i = 0; i < point.length; i++)
                centroid[i] += point[i];
        }

        //Then we average them
        for (int i = 0; i < centroid.length; i++)
            centroid[i] = centroid[i] / points.length;
        return centroid;
    }

    /**
     * Computes the cross-product of two 3D vectors
     * @param v  the first vector
     * @param w  the second vector
     * @return the new vector, containing the result of the cross product
     * @exception IllegalArgumentException We can only compute the cross product in 3D
     * @exception IllegalArgumentException The vectors do not have the same dimensions
     */
    public static double[] crossProduct3(double[] v, double[] w){
        try {
            if (v.length != 3) throw new IllegalArgumentException("We can only compute the cross product in 3D");
            if(v.length != w.length )
                throw new IllegalArgumentException("The vectors do not have the same dimensions");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        return new double[]{
                v[1] * w[2] - v[2] * w[1],
                v[2] * w[0] - v[0] * w[2],
                v[0] * w[1] - v[1] * w[0]
        };
    }

    /**
     * Computes the magnitude (length) of the vector
     * @param v  the vector
     * @return the magnitude of v
     */
    public static double magnitude(double[] v) {
        double sqSum = 0;
        for(int i = 0; i < v.length; i++)  sqSum += v[i] * v[i];
        return Math.sqrt(sqSum);
    }

    /**
     * Computes the unit vector of v
     * @param v  the vector
     * @return the unit vector of v
     */
    public static double[] normalize(double[] v){
        double l = magnitude(v);
        double[] w = new double[v.length];
        for(int i = 0; i < v.length; i++)    w[i] = v[i]/l;
        return w;
    }

    /**
     * Computes the dot product of two vectors
     * @param v  the first vector
     * @param w  the second vector
     * @return   the dot product of v and u
     * @exception IllegalArgumentException The vectors do not have the same dimensions
     */
    public static double dotProduct(double[] v, double[] w){
        try {
            if(v.length != w.length )
                throw new IllegalArgumentException("The vectors do not have the same dimensions");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        double sum = 0;
        for(int i = 0; i < v.length; i++)   sum += v[i] * w[i];
        return sum;
    }

    /**
     * Computes the subtraction of two vectors
     * @param v  the first vector
     * @param w  the second vector
     * @return   the difference between v and w
     * @exception IllegalArgumentException The vectors do not have the same dimensions
     */
    public static double[] vectorSubtract(double[] v, double[] w){
        try {
            if(v.length != w.length )
                throw new IllegalArgumentException("The vectors do not have the same dimensions");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        double[] x = new double[v.length];
        for(int i = 0; i < v.length; i++)   x[i] = v[i] - w[i];
        return x;
    }

    /**
     * Computes the sum of two vectors
     * @param v  the first vector
     * @param w  the second vector
     * @return   the sum of v and w
     * @exception IllegalArgumentException The vectors do not have the same dimensions
     */
    public static double[] vectorSum(double[] v, double[] w){
        try {
            if(v.length != w.length )
                throw new IllegalArgumentException("The vectors do not have the same dimensions");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        double[] x = new double[v.length];
        for(int i = 0; i < v.length; i++)   x[i] = v[i] + w[i];
        return x;
    }


    /**
     * Computes the euclidean distance between two points in euclidean space.
     * @param a  the first point (a vector of euclidean coordinates).
     * @param b  the second point (^^).
     * @return the euclidean distance between the two points.
     * @exception IllegalArgumentException The Cartesian coordinates vectors do not have the same dimensions.
     */
    public static double euclideanDistance( double[] a,  double[] b){
        try {
            if(a.length != b.length )
                throw new IllegalArgumentException("The Cartesian coordinates vectors do not have the same dimensions");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        double squared = 0.0;
        for (int i = 0; i < a.length; i++)
            squared += (b[i] - a[i]) * (b[i] - a[i]);
        return Math.sqrt(squared);
    }

    /**
     * Computes the Manhattan distance between two points
     * @param a  the first point (a vector of euclidean coordinates).
     * @param b  the second point (^^).
     * @return the euclidean distance between the two points.
     * @exception IllegalArgumentException The Cartesian coordinates vectors do not have the same dimensions.
     */
    public static double manhattanDistance( double[] a,  double[] b){
        try {
            if(a.length != b.length )
                throw new IllegalArgumentException("The Cartesian coordinates vectors do not have the same dimensions");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        double total = 0;
        for (int i = 0; i < a.length; i++)
            total+= Math.abs(a[i] - b[i]);
        return total;
    }

    /**
     * Computes a rotation matrix R to rotate the vector ab to the vector ac
     * BASED ON: https://en.wikipedia.org/wiki/Rotation_matrix#Vector_to_vector_formulation
     * @param a the origin point
     * @param b the first target point
     * @param c the second target point
     * @return the rotation matrix R such that R*ab = ac
     */
    public static Matrix rotationMatrix(double[] a, double[] b, double[] c){
        //get ab and ac as unit vectors
        double[] ab = normalize(vectorSubtract(b,a));
        double[] ac = normalize(vectorSubtract(c,a));

        //column vectors as matrix
        Matrix AB = new Matrix(ab,3);
        Matrix AC = new Matrix(ac,3);

        //row vectors as matrix
        Matrix ABT = new Matrix(ab,1);
        Matrix ACT = new Matrix(ac,1);

        //compute yx and xy
        Matrix yx = AC.times(ABT);
        Matrix xy = AB.times(ACT);
        Matrix yxMxy = yx.minus(xy);
        Matrix yxMxy2 = yxMxy.times(yxMxy);

        double mul = 1 / (1 + dotProduct(ab, ac));

        Matrix last = yxMxy2.times(mul);

        //put everything together and compute R

        return Matrix.identity(3,3).plus(yxMxy).plus(last);
    }

    /**
     * Converts a ("pure") 3D rotation matrix into a quaternion
     * @param r the rotation matrix
     * @return The quaternion
     * @throws IllegalArgumentException The {@link Matrix} r does not have the right dimensions
     * @throws IllegalArgumentException The {@link Matrix} r is not a "pure" rotation matrix
     * @throws IllegalArgumentException The {@link Matrix} r is not special orthogonal
     */
    public static double[] rotationToQuaternion(Matrix r){
        try {
            double roundingError = 0.0001;
            if(r.getColumnDimension() != 3 || r.getRowDimension() != 3)
                throw new Exception("The Matrix r does not have the right dimensions");
            if(!r.times(r.transpose()).epsEqual(Matrix.identity(3,3), roundingError))
                throw new Exception("The Matrix r is not a \"pure\" rotation matrix");

            if(!Maths.epsEqual(r.determinant(), 1.0, roundingError))
                throw new Exception("The Matrix r is not special orthogonal");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        double tr = r.trace();
        double qw,qx,qy,qz;

        double m00 = r.get(0,0);
        double m01 = r.get(0,1);
        double m02 = r.get(0,2);
        double m10 = r.get(1,0);
        double m11 = r.get(1,1);
        double m12 = r.get(1,2);
        double m20 = r.get(2,0);
        double m21 = r.get(2,1);
        double m22 = r.get(2,2);

        if (tr > 0) {
            double S = Math.sqrt(tr+1.0) * 2; // S=4*qw
            qw = 0.25 * S;
            qx = (m21 - m12) / S;
            qy = (m02 - m20) / S;
            qz = (m10 - m01) / S;
        } else if ((m00 > m11)&(m00 > m22)) {
            double S = Math.sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx
            qw = (m21 - m12) / S;
            qx = 0.25 * S;
            qy = (m01 + m10) / S;
            qz = (m02 + m20) / S;
        } else if (m11 > m22) {
            double S = Math.sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
            qw = (m02 - m20) / S;
            qx = (m01 + m10) / S;
            qy = 0.25 * S;
            qz = (m12 + m21) / S;
        } else {
            double S = Math.sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
            qw = (m10 - m01) / S;
            qx = (m02 + m20) / S;
            qy = (m12 + m21) / S;
            qz = 0.25 * S;
        }

        return new double[]{qw,qx,qy,qz};
    }

    /**
     * Rotate a vector around another vector in 3D
     * @param v            the vector to rotate
     * @param u            the axis to rotate around (unit vector)
     * @param angle        the rotation angle (in radians)
     * @exception IllegalArgumentException We can only do the vector rotation in 3D
     * @exception IllegalArgumentException The vectors do not have the same dimensions
     */
    public static double[] rodriguesGibbs(double[] v, double[] u, double angle){
        try {
            if (v.length != 3) throw new IllegalArgumentException("We can only do the vector rotation in 3D");
            if(v.length != u.length )
                throw new IllegalArgumentException("The vectors do not have the same dimensions");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        double[] uxv = crossProduct3(u, v);
        double udv = dotProduct(u,v);

        double[] w = new double[3];
        for(int i = 0; i < 3; i++){
            w[i] = v[i] * Math.cos(angle) + uxv[i] * Math.sin(angle) + udv * u[i] * (1 - Math.cos(angle));
        }

        return w;
    }

    /**
     * Converts between cosine and sine (should work both ways)
     * Note that we lose the sign of the value...
     * @param other  the cosine or sine value of an angle
     * @return        the converted value
     */
    public static double convertCosSin(double other){
        return Math.sqrt(1 - other * other);
    }

    /**
     * Computes the cosine angle abc between three points A, B, C using the cosine rule
     * @param  ab          the distance between A and B
     * @param  bc          the distance between B and C
     * @param  ac          the distance between A and C (this side is opposite the angle abc)
     * @return             the cosine of the angle abc (the angle at point B)
     */
    public static double cosSSS(double ab, double bc, double ac){
        return (ab*ab + bc * bc - ac * ac)  / (2 * ab * bc);
    }

    /**
     * Computes the angle abc between three points A, B, C using the cosine rule
     * @param  ab          the distance between A and B
     * @param  bc          the distance between B and C
     * @param  ac          the distance between A and C (this side is opposite the angle abc)
     * @return             the angle abc (the angle at point B) in radians
     */
    public static double solveSSS(double ab, double bc, double ac)
    {
        return Math.acos(cosSSS(ab,bc,ac));
    }

    /**
     * Computes the distance ac, given two sides and the angle between them (SAS)
     *
     * @param  ab          the distance between A and B
     * @param  bc          the distance between B and C
     * @param  abc         the angle abc (the angle at point B) in radians
     * @return             the distance ac
     */
    public static double solveSAS(double ab, double bc, double abc)
    {
        return Math.sqrt(ab * ab + bc * bc - 2 * ab * bc * Math.cos(abc));
    }

    /**
     * Computes the angle between a line and the x-axis. The line is defined by two points, a and b.
     * @param a the first point on the line
     * @param b the second point on the line
     * @return the angle (in radians)
     */
    public static double xAngle(double[] a, double[] b){
        try {
            if (a.length != 2 || b.length != 2) throw new IllegalArgumentException("This method only works in 2D");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        double angle = Math.atan2(b[1]-a[1], b[0] - a[0]);
        if(angle < 0) return xAngle(b,a);
        else return angle;
    }


    /**
     * Computes the unit normal vector given 3 points in a plane
     *
     * @param  a          the coordinates of point a
     * @param  b          the coordinates of point b
     * @param  c          the coordinates of point c
     * @return            the normal vector of the plane defined by a, b and c
     * @exception IllegalArgumentException We can only compute the plane in 3D
     * @exception IllegalArgumentException The three points do not have the same dimensions
     */
    public static double[] planeNormal(double[] a, double[] b, double[] c){
        try {
            if (a.length != 3) throw new IllegalArgumentException("We can only compute the plane in 3D");
            if(a.length != b.length || b.length != c.length)
                throw new IllegalArgumentException("The three points do not have the same dimensions");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        int dim = 3;

        //vectors AB, BC and unit vector BC
        double[] AB = new double[dim];
        double[] AC = new double[dim];
        for(int i = 0; i < dim; i++) {
            AB[i] = b[i] - a[i];
            AC[i] = c[i] - a[i];
        }

        double[] n = crossProduct3(AB, AC);
        return normalize(n);
    }

    /**
     * Given three (consecutive) distances between four points, as well as two angles and a torsion angle,
     * computes the distance between the first and the last point of the sequence.
     * Using: https://math.stackexchange.com/questions/4175682/distance-between-two-points-in-two-different-planes-given-the-dihedral-angle/4181249
     * This method is based on the first answer on stackoverflow
     * @param ab             the distance between a and b
     * @param bc             the distance between b and c
     * @param cd             the distance between c and d
     * @param abc            the angle between ab and bc (in radians)
     * @param bcd            the angle between bc and cd (in radians)
     * @param dihedral       the measured dihedral angle (in radians)
     * @return               the distance between a and d based on the dihedral angle
     */
    public static double dihedralDistance(double ab, double bc, double cd, double abc, double bcd, double dihedral)
    {
        //To calculate the distance we do the following:
        //We set the position of b to (0,0,0)
        //We compute the position of a using sin/cos
        //We compute  the position of d similarly, assuming that the torsion angle is 0
        //Finally, we rotate d around the x-axis using a rotation matrix.
        Matrix a = new Matrix(3,1);
        a.add(0,0, ab * Math.cos(abc));  a.add(1,0, ab * Math.sin(abc));  a.add(2,0, 0);

        Matrix d = new Matrix(3,1);
        d.add(0,0, -cd * Math.cos(bcd) + bc);  d.add(1,0, cd * Math.sin(bcd));  d.add(2,0, 0);

        if(Math.toDegrees(dihedral) != 0){
            Matrix rotation = new Matrix(3,3);
            rotation.add(0,0,1);
            rotation.add(0,1,0);
            rotation.add(0,2,0);
            rotation.add(1,0,0);
            rotation.add(2,0,0);
            rotation.add(1,1, Math.cos(dihedral));
            rotation.add(1,2, Math.sin(dihedral));
            rotation.add(2,1,-Math.sin(dihedral));
            rotation.add(2,2, Math.cos(dihedral));

            d = rotation.times(d);
        }

        return euclideanDistance(a.getColumn(0), d.getColumn(0));
    }

    /**
     * Given three (consecutive) distances between four points, as well as two angles and a torsion angle,
     * computes the distance between the first and the last point of the sequence.
     * Using: https://math.stackexchange.com/questions/4175682/distance-between-two-points-in-two-different-planes-given-the-dihedral-angle/4181249
     * This method is based on the second answer on stackoverflow
     * @param ab             the distance between a and b
     * @param bc             the distance between b and c
     * @param cd             the distance between c and d
     * @param abc            the angle between ab and bc (in radians)
     * @param bcd            the angle between bc and cd (in radians)
     * @param dihedral       the measured dihedral angle (in radians)
     * @return               the distance between a and d based on the dihedral angle
     */
    public static double dihedralDistance2(double ab, double bc, double cd, double abc, double bcd, double dihedral){
        dihedral = Math.abs(dihedral);

        double abcSup = Math.toRadians(180 - Math.toDegrees(abc));
        double bcdSup = Math.toRadians(180 - Math.toDegrees(bcd));

        double[] a = new double[]{
                ab * Math.sin(abcSup) * Math.cos(dihedral),
                0,
                ab * Math.sin(abcSup) * Math.sin(dihedral)};

        double[] d = new double[]{
                cd * Math.sin(bcdSup),
                ab * Math.cos(abcSup) + bc + cd * Math.cos(bcdSup),
                0};
        return euclideanDistance(a, d);
    }


    /**
     * Computes the (SIGNED) torsion angle between the planes of a and d, given coordinates of all four points
     * Using: https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
     * @param a             the coordinates of points a
     * @param b             the coordinates of points b
     * @param c             the coordinates of points c
     * @param d             the coordinates of points d
     * @exception IllegalArgumentException We can only find the dihedral angle between objects in 3D space
     * @exception IllegalArgumentException The four points do not have the same dimensions
     * @return              the dihedral angle in radians
     */
    public static double signedDihedralAngle(double[] a, double[] b, double[] c, double[] d){
        try {
            if (a.length != 3) throw new IllegalArgumentException("We can only find the dihedral angle between objects in 3D space");
            if(a.length != b.length || b.length != c.length || c.length != d.length)
                throw new IllegalArgumentException("The four points do not have the same dimensions");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        double[] n1 = planeNormal(a, b, c);
        double[] n2 = planeNormal(b,c,d);
        double[] BC = vectorSubtract(c,b);
        double[] m1 = crossProduct3(n1, normalize(BC));
        double x = dotProduct(n1,n2);
        double y = dotProduct(m1,n2);
        return -Math.atan2(y,x);
    }

    /**
     * Computes the cosine of the dihedral angle created by the planes abc and bcd.
     * @param ab the distance from a to b
     * @param bc the distance from b to c
     * @param ac the distance from a to c
     * @param ad the distance from a to d
     * @param bd the distance from b to d
     * @param cd the distance from c to d
     * @return the cosine of the omega value if valid, -1 otherwise
     */
    public static double cosOmega(double ab, double bc, double ac, double ad, double bd, double cd) {
        //Use the cosine law to compute several angle values
        double cosABD = cosSSS(ab, bd, ad); //corresponds to "a" in mdjeep cosomega
        double cosCBD = cosSSS(bc, bd, cd); //corresponds to "b" in mdjeep cosomega
        double cosABC = cosSSS(ab, bc, ac); //corresponds to "c" in mdjeep cosomega
        double sinCBD = convertCosSin(cosCBD); //corresponds to "e" in mdjeep cosomega
        double sinABC = convertCosSin(cosABC); //corresponds to "f" in mdjeep cosomega
        double cosOmega = (cosABD - cosCBD * cosABC) / (sinCBD * sinABC);

        double tolerance = 0.1;
        if(cosOmega > 1 && cosOmega < 1 + tolerance){
            return 1;
        }
        if(cosOmega < -1 && cosOmega > -1 - tolerance){
            return -1;
        }
        return cosOmega;
    }


    /**
     * Computes the torsion angle between the planes of a and d, using pairwise distances
     * USING: The discretizable molecular distance geometry problem, C.Lavor, L.Liberti, N. Maculan and A.M.
     * @param ab            the distance from a to b
     * @param bc            the distance from b to c
     * @param ac            the distance from a to c
     * @param ad            the distance from a to d
     * @param bd            the distance from b to d
     * @param cd            the distance from c to d
     * @return              the dihedral angle of the planes abc and bcd
     */
    public static double dihedralAngle(double ab, double bc, double ac, double ad, double bd, double cd) {
        double cosOmega = cosOmega(ab, bc, ac, ad, bd, cd);
        double sinOmega = convertCosSin(cosOmega);
        return Math.atan2(sinOmega, cosOmega);
    }

    /**
     * Computes the torsion angle between the planes of a and d, using pairwise distances
     * Performs a stepwise search in the interval [adL, adU] in order to find a feasible dihedral angle
     * USING: The discretizable molecular distance geometry problem, C.Lavor, L.Liberti, N. Maculan and A.M.
     * @param ab            the distance from a to b
     * @param bc            the distance from b to c
     * @param ac            the distance from a to c
     * @param adL           the lowerbound on distance from a to d
     * @param adL           the upperbound on distance from a to d
     * @param bd            the distance from b to d
     * @param cd            the distance from c to d
     * @param mode          the point to start from, 0 = lowerbound, 1 = upperbound
     * @param step          the stepsize (eps in C-MDJeep)
     * @exception IllegalArgumentException Mode should be either 0 or 1
     * @return              the lowerbound or upperbound on the (cosine of) the dihedral angle of the planes abc and bcd given an interval distance ad, or null when there is no feasible angle
     */
    public static Double feasibleDihedralAngle(double ab, double bc, double ac, double adL, double adU, double bd, double cd, int mode, double step){
        try{
            if(mode != 0 && mode != 1)
                throw new IllegalArgumentException("Mode should be either 0 or 1");
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        double r = mode; //r starts at 0.0 or 1.0
        double cosOmega = -2; //initialize outside [-1,1] interval

        // computing cosine of first feasible angle
        while (r >= 0.0 && r <= 1.0 && (cosOmega < -1.0 || cosOmega > 1.0))
        {
            double ad = adL + r * (adU - adL);
            cosOmega = cosOmega(ab, bc, ac, ad, bd, cd);

            // preparing for next loop?
            if (mode == 0)  r += step;
            if (mode == 1)  r -= step;
        };



        // if feasible angle was not found, return null
        if (cosOmega < -1.0 || cosOmega > 1.0) return null;
        double sinOmega = convertCosSin(cosOmega);
        return Math.atan2(sinOmega, cosOmega);
    }

    /**
     * Computes the coordinates of the fourth point d point given the coordinates of the other 3 points (a, b, c)
     * as well as the distances bd and cd and the torsion angle
     * Based on the Rodrigues-Gibbs Formulation (RG) of the Two-Step vector rotation method
     * https://pubmed.ncbi.nlm.nih.gov/15898109/
     * @param a             the coordinates of points a
     * @param b             the coordinates of points b
     * @param c             the coordinates of points c
     * @param bc            the distance between b and c, can be null (then we simply compute it
     * @param bd            the distance between b and d
     * @param cd            the distance between c and d
     * @param dihedral      the measured (signed) dihedral angle (in radians)
     * @exception IllegalArgumentException We can only use the dihedral angle between objects in 3D space
     * @exception IllegalArgumentException The three points do not have the same dimensions
     * @return              the coordinates of d
     */
    public static double[] dihedralCoordinates(double[] a,double[] b,double[] c, Double bc, double bd,double cd, double dihedral) {
        try {
            if (a.length != 3)
                throw new IllegalArgumentException("We can only use the dihedral angle between objects in 3D space");
            if (a.length != b.length || b.length != c.length)
                throw new IllegalArgumentException("The three points do not have the same dimensions");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        int dim = 3;

        //STEP 0: initalize values
        bc = bc == null ? euclideanDistance(b, c) : bc;

        //compute angle bcd
        double bcd = solveSSS(bc, cd, bd);
        bcd = Math.toRadians(180) - bcd;

        //vector BC
        double[] uBC = normalize(vectorSubtract(c,b));

        //STEP 1: compute coordinates of d0, which is cd away from point c extending along the BC axis
        //For simplicity, we set C as the origin
        double[] d0 = new double[dim];
        for (int i = 0; i < dim; i++) d0[i] = cd * uBC[i];

        //normal vector
        double[] nu = planeNormal(b, c, a);

        //STEP 2: compute d1, we rotate d0 around the normal by the angle bcd
        double[] d1 = rodriguesGibbs(d0, nu, bcd);

        //STEP 3: compute d2, we rotate d1 around uBC by the torsion angle
        double[] d2 = rodriguesGibbs(d1, uBC, dihedral);
        return vectorSum(d2, c);
    }

    /**
     * Computes the coordinates of the fourth point d point given the coordinates of the other 3 points (a, b, c)
     * as well as the distances bd and cd and the torsion angle
     * Based on the Rodrigues-Gibbs Formulation (RG) of the Two-Step vector rotation method
     * https://pubmed.ncbi.nlm.nih.gov/15898109/
     * @param a             the coordinates of points a
     * @param b             the coordinates of points b
     * @param c             the coordinates of points c
     * @param bcd           the vector angle bcd in radians
     * @param cd            the distance between c and d
     * @param dihedral      the measured (signed) dihedral angle (in radians)
     * @exception IllegalArgumentException We can only use the dihedral angle between objects in 3D space
     * @exception IllegalArgumentException The three points do not have the same dimensions
     * @return              the coordinates of d
     */
    public static double[] dihedralCoordinates(double[] a,double[] b,double[] c, double cd, double bcd, double dihedral) {
        try {
            if (a.length != 3)
                throw new IllegalArgumentException("We can only use the dihedral angle between objects in 3D space");
            if (a.length != b.length || b.length != c.length)
                throw new IllegalArgumentException("The three points do not have the same dimensions");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        int dim = 3;
        bcd = Math.toRadians(180) - bcd;

        //vector BC
        double[] uBC = normalize(vectorSubtract(c,b));

        //normal vector
        double[] nu = planeNormal(b, c, a);

        //STEP 1: compute coordinates of d0, which is cd away from point c extending along the BC axis
        //For simplicity, we set C as the origin
        double[] d0 = new double[dim];
        for (int i = 0; i < dim; i++) d0[i] = cd * uBC[i];

        //STEP 2: compute d1, we rotate d0 around the normal by the angle bcd
        double[] d1 = rodriguesGibbs(d0, nu, bcd);

        //STEP 3: compute d2, we rotate d1 around uBC by the torsion angle
        double[] d2 = rodriguesGibbs(d1, uBC, dihedral);
        return vectorSum(d2, c);
    }

    /**
     * Computes the coordinates of the fourth point d point given the coordinates of the other 3 points (a, b, c)
     * as well as the distances bd and cd and the cosine of the torsion angle omega
     * based on:
     * D.S. Gonçalves, A. Mucherino, Discretization Orders and Efficient Computation of Cartesian Coordinates for Distance Geometry, Optimization Letters 8(7), 2111-2125, 2014.
     * @param a             the coordinates of points a
     * @param b             the coordinates of points b
     * @param c             the coordinates of points c
     * @param bc            the distance between b and c, can be null (then we simply compute it
     * @param bd            the distance between b and d
     * @param cd            the distance between c and d
     * @param cosOmega      the cosine of the dihedral angle (omega)
     * @param sinOmega      the sine of the dihedral angle (omega) NOTE: the sign of sinOmega decides which of the two intersections we are choosing!
     * @exception IllegalArgumentException We can only use the dihedral angle between objects in 3D space
     * @exception IllegalArgumentException The three points do not have the same dimensions
     * @return              the coordinates of d
     */
    public static double[] dihedralCoordinates(double[] a,double[] b,double[] c, double bc, double bd,double cd, double cosOmega, double sinOmega) {
        try {
            if (a.length != 3)
                throw new IllegalArgumentException("We can only use the dihedral angle between objects in 3D space");
            if (a.length != b.length || b.length != c.length)
                throw new IllegalArgumentException("The three points do not have the same dimensions");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        double[][] uMatrix = uMatrix(a,b,c);
        double cosTheta = cosSSS(bc,cd,bd);
        double sinTheta = convertCosSin(cosTheta);

        return dihedralCoordinates(uMatrix, c, cd, cosOmega, sinOmega, cosTheta, sinTheta);

    }

    public static double[] dihedralCoordinates(double[][] uMatrix, double[] c,double cd, double cosOmega, double sinOmega, double cosTheta, double sinTheta) {

        // computing vector a (depends on angles)
        double a0 = -cd*cosTheta;
        double a1 =  cd*sinTheta*cosOmega;
        double a2 =  cd*sinTheta*sinOmega;

        // generation of the coordinates
        double[] d = new double[3];
        d[0] = c[0] + a0 * uMatrix[0][0] + a1 * uMatrix[1][0] + a2 * uMatrix[2][0];
        d[1] = c[1] + a0 * uMatrix[0][1] + a1 * uMatrix[1][1] + a2 * uMatrix[2][1];
        d[2] = c[2] + a0 * uMatrix[0][2] + a1 * uMatrix[1][2] + a2 * uMatrix[2][2];
        return d;
    }


        public static double[][] uMatrix(double[] a, double[] b, double[] c){
        //compute the U matrix
        double[] v1 = vectorSubtract(c,b); //the vector from b to c
        double[] v2 = vectorSubtract(a,b); //the vector from b to a

        double[][] uMatrix = new double[3][];

        // x axis (first column)
        uMatrix[0] = normalize(v1);

        //z axis (third column)
        uMatrix[2] = normalize(crossProduct3(v1,v2));

        //y axis (second column)
        uMatrix[1] = normalize(crossProduct3(uMatrix[2],uMatrix[0]));

        return uMatrix;
    }

    /**
     * Translates A and B and rotates B such that A and B are as best aligned as possible
     * https://en.wikipedia.org/wiki/Kabsch_algorithm
     *
     * @param  A   the first set of points
     * @param  B   the second set of points
     * @exception IllegalArgumentException The point sets must be of the same length.
     * @exception IllegalArgumentException The point sets must contain at least 3 points (A and B must have at least 3 columns).
     * @exception IllegalArgumentException The points must me in dimension 3 (A and B must contain 3 rows).
     */
    public static void kabschAlignment(Matrix A, Matrix B){
        try{
            if(A.getRowDimension() != B.getRowDimension())
                throw new IllegalArgumentException("The point sets must be of the same length.");
//            if(A.getRowDimension() < 3)
//                throw new IllegalArgumentException("The point sets must contain at least 3 points (A and B must have at least 3 columns).");
//            if(A.getColumnDimension() != 3 || B.getColumnDimension() != 3 )
//                throw new IllegalArgumentException("The points must me in dimension 3 (A and B must contain 3 columns).");
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        //Set the centroid of both point sets to the origin
        double[] centroidA = centroid(A.getArray());
        double[] centroidB = centroid(B.getArray());

        //Set the centroids of both A and B to 0,0,0 (subtract the centroid from each of the points)
        for(int i = 0; i < A.getRowDimension(); i++){
            for(int j = 0; j < A.getColumnDimension(); j++)
            {
                A.set(i,j, A.get(i,j) - centroidA[j]);
                B.set(i,j, B.get(i,j) - centroidB[j]);
            }
        }



        //The covariance matrix H
        Matrix H = A.transpose().times(B);

        Matrix U = H.SVD_getU();
        Matrix V = H.SVD_getV();

        //Reflection exception
        boolean reflect = U.times(V).determinant() < 0; //U.determinant() * V.determinant() == -1
        if(reflect)
            U.times(-1);

        Matrix rotation = V.times(U.transpose());

        //Rotate B
        Matrix BRotated = B.times(rotation);
        for(int i = 0; i < B.getRowDimension(); i++)
            for(int j = 0; j < B.getColumnDimension(); j++)
                B.set(i,j, BRotated.get(i,j));
    }

    /**
     * Computes the RMSD between two point sets in 3D
     *
     * @param  A   the first set of points, represented by a matrix
     * @param  B   the second set of points, represented by a matrix
     * @return     the root mean square deviation
     */
    public static double RMSD3D(Matrix A, Matrix B){
        double dev = 0;
        //Sum euclidean distances and divide by n
        for(int i = 0; i < A.getRowDimension(); i++) {
            dev += Math.pow(euclideanDistance(A.getRow(i), B.getRow(i)), 2);
        }

        //Take the square root
        return Math.sqrt(dev/ A.getRowDimension());
    }
}
