/*
 * Useful XYZ processing tool
 */
package cmtlab;

import Jama.Matrix;
import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;

/**
 *
 * @author Sherif Abdulkader Tawfik
 * sherif.tawfic@gmail.com
 * sherif.abbas@rmit.edu.au
 */
public class XYZ {

    public String systemName;
    public String atom[];
    public double x[];
    public double y[];
    public double z[];
    public int N;

    public int numSpecies;
    public int speciesIndex[];
    public int speciesMassNumber[];
    public double speciesAtomicMass[];
    public String speciesLabel[];
    public ArrayList v_speciesLabel = new ArrayList();

    public double doubleTag1;

    public XYZ(int N) {
        this.N = N;
        atom = new String[N];
        x = new double[N];
        y = new double[N];
        z = new double[N];
    }

    public XYZ(String atom, double x, double y, double z) {
        this(1);
        this.atom[0] = new String(atom);
        this.x[0] = x;
        this.y[0] = y;
        this.z[0] = z;

        setup();

    }

    public XYZ(XYZ xyz) {
        this(xyz.N);
        this.atom = xyz.atom.clone();
        this.x = xyz.x.clone();
        this.y = xyz.y.clone();
        this.z = xyz.z.clone();
        if (xyz.systemName != null) {
            this.systemName = new String(xyz.systemName);
        }
        this.doubleTag1 = xyz.doubleTag1;
        setup();
    }

    public XYZ(String fileName) {
        BufferedReader in = null;
        try {

            in = new BufferedReader(new FileReader(new File(fileName)));

            String line = in.readLine();

            N = Integer.valueOf(line.replace(".", "").trim());
            atom = new String[N];
            x = new double[N];
            y = new double[N];
            z = new double[N];

            in.readLine();

            for (int i = 0; i < N; i++) {
                line = in.readLine().trim();
                String lineCont[] = line.split("\\s+");

                atom[i] = lineCont[0];

                x[i] = Double.valueOf(lineCont[1]);
                y[i] = Double.valueOf(lineCont[2]);
                z[i] = Double.valueOf(lineCont[3]);
                setup();
            }
        } catch (IOException e) {
            e.printStackTrace();

        } finally {
            if (in != null) {
                try {
                    in.close();
                } catch (IOException ioe) {
                    ioe.printStackTrace();
                }
            }
        }
    }

    public void setup() {

        v_speciesLabel = new ArrayList();
        for (int i = 0; i < N; i++) {
            if (!v_speciesLabel.contains(atom[i])) {
                v_speciesLabel.add(atom[i]);
            }
        }

        numSpecies = v_speciesLabel.size();
        speciesLabel = new String[numSpecies];

        for (int i = 0; i < numSpecies; i++) {
            speciesLabel[i] = (String) v_speciesLabel.get(i);
        }

        speciesIndex = new int[numSpecies];
        speciesMassNumber = new int[numSpecies];
        speciesAtomicMass = new double[numSpecies];
        for (int i = 0; i < numSpecies; i++) {
            speciesIndex[i] = i + 1;
            speciesMassNumber[i] = AtomicData.getMassNumber(speciesLabel[i]);
            speciesAtomicMass[i] = AtomicData.getAtomicMass(speciesLabel[i]);
        }
    }

    public int getSpeciesListIndex(String label) {

        for (int i = 0; i < numSpecies; i++) {

            if (speciesLabel[i].equals(label)) {
                return speciesIndex[i];
            }
        }
        System.out.println("Cannot find atom " + label);
        System.exit(1);
        return -1;
    }

    public String getSystemName() {
        return systemName;
    }

    public int howManyOf(String a) {
        int count = 0;
        for (int i = 0; i < N; i++) {
            if (atom[i].equals(a)) {
                count++;
            }
        }
        return count;
    }

    public double distance(int a, int b) {
        a = a - 1;
        b = b - 1;
        return Math.sqrt((x[a] - x[b]) * (x[a] - x[b]) + (y[a] - y[b]) * (y[a] - y[b]) + (z[a] - z[b]) * (z[a] - z[b]));
    }

    public double distance2D(int a, int b) {
        a = a - 1;
        b = b - 1;
        return Math.sqrt((x[a] - x[b]) * (x[a] - x[b]) + (y[a] - y[b]) * (y[a] - y[b]));
    }

    public double distance(int a, double b[]) {
        a = a - 1;
        return Math.sqrt((x[a] - b[0]) * (x[a] - b[0]) + (y[a] - b[1]) * (y[a] - b[1]) + (z[a] - b[2]) * (z[a] - b[2]));
    }

    public double distance2D(int a, double b[]) {
        a = a - 1;
        return Math.sqrt((x[a] - b[0]) * (x[a] - b[0]) + (y[a] - b[1]) * (y[a] - b[1]));
    }

    public double distance(int[] a, int[] b) {
        double min = 10, dist = 0;
        for (int k = 0; k < a.length; k++) {
            for (int i = 0; i < b.length; i++) {
                dist = Math.sqrt((x[a[k] - 1] - x[b[i] - 1]) * (x[a[k] - 1] - x[b[i] - 1]) + (y[a[k] - 1] - y[b[i] - 1]) * (y[a[k] - 1] - y[b[i] - 1]) + (z[a[k] - 1] - z[b[i] - 1]) * (z[a[k] - 1] - z[b[i] - 1]));
                if (min > dist) {
                    min = dist;
                }
            }
        }
        return min;
    }

    public double distanceAxis(int axis, int[] a, int[] b) {
        double min = 10000, dist = 0;
        if (axis == 3) {
            for (int k = 0; k < a.length; k++) {
                for (int i = 0; i < b.length; i++) {
                    dist = Math.abs(z[a[k] - 1] - z[b[i] - 1]);
                    if (min > dist) {
                        min = dist;
                    }
                }
            }
        }
        return min;
    }

    public double distance(int a, int[] b) {
        a = a - 1;
        double min = 10, dist = 0;
        for (int i = 0; i < b.length; i++) {
            dist = Math.sqrt((x[a] - x[b[i] - 1]) * (x[a] - x[b[i] - 1]) + (y[a] - y[b[i] - 1]) * (y[a] - y[b[i] - 1]) + (z[a] - z[b[i] - 1]) * (z[a] - z[b[i] - 1]));
            if (min > dist) {
                min = dist;
            }
        }
        return min;
    }

    public double distance(int a, XYZ b) {
        a = a - 1;
        double min = 10, dist = 0;
        for (int i = 1; i <= b.N; i++) {
            dist = Math.sqrt((x[a] - b.x[i - 1]) * (x[a] - b.x[i - 1]) + (y[a] - b.y[i - 1]) * (y[a] - b.y[i - 1]) + (z[a] - b.z[i - 1]) * (z[a] - b.z[i - 1]));
            if (min > dist) {
                min = dist;
            }
        }
        return min;
    }

    public XYZ replace(int index, String newLabel) {
        atom[index] = new String(newLabel);
        return this;
    }

    public int getSpeciesIndex(String label) {

        for (int i = 0; i < N; i++) {

            if (atom[i].equals(label)) {
                return i;
            }
        }
        System.out.println("Cannot find atom " + label);
        System.exit(1);
        return -1;
    }

    public void writeTo(String fileName) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
            NumberFormat formatter = new DecimalFormat("#0.00000000000000000000");

            out.write("" + N + "\n\n");
            for (int i = 0; i < N; i++) {
                out.write(atom[i] + "\t" + formatter.format(x[i]) + "\t" + formatter.format(y[i]) + "\t" + formatter.format(z[i]) + "\n");
                out.flush();
            }
        } catch (IOException e) {
            e.printStackTrace();

        }
    }

    public void writeToGAUSSIANWithVariables(String fileName) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
            NumberFormat formatter = new DecimalFormat("#0.00000000000000000000");

            for (int i = 0; i < N; i++) {
                out.write(atom[i] + "\tA" + i + "1" + "\tA" + i + "2" + "\tA" + i + "3" + "\n");
                out.flush();
            }
            for (int i = 0; i < N; i++) {
                out.write("A" + i + "1\t" + formatter.format(x[i]) + "\n");
                out.write("A" + i + "2\t" + formatter.format(y[i]) + "\n");
                out.write("A" + i + "3\t" + formatter.format(z[i]) + "\n");
                out.flush();
            }
        } catch (IOException e) {
            e.printStackTrace();

        }
    }

    public XYZ centerX() {

        double minX = Utilities.min(x);
        double maxX = Utilities.max(x);
        double avX = (minX + maxX) / 2;
        for (int i = 0; i < N; i++) {
            x[i] = x[i] - avX;

        }
        return this;
    }

    public XYZ centerY() {
        double miny = Utilities.min(y);
        double maxy = Utilities.max(y);
        double avy = (miny + maxy) / 2;
        for (int i = 0; i < N; i++) {
            y[i] = y[i] - avy;

        }
        return this;
    }

    public XYZ centerZ() {
        double minz = Utilities.min(z);
        double maxz = Utilities.max(z);
        double avz = (minz + maxz) / 2;
        for (int i = 0; i < N; i++) {
            z[i] = z[i] - avz;

        }
        return this;
    }

    public XYZ center() {
        double minX = Utilities.min(x);
        double maxX = Utilities.max(x);
        double minY = Utilities.min(y);
        double maxY = Utilities.max(y);
        double minZ = Utilities.min(z);
        double maxZ = Utilities.max(z);
        double avX = (minX + maxX) / 2;
        double avY = (minY + maxY) / 2;
        double avZ = (minZ + maxZ) / 2;

        for (int i = 0; i < N; i++) {
            x[i] = x[i] - avX;
            y[i] = y[i] - avY;
            z[i] = z[i] - avZ;
        }
        return this;
    }

    public XYZ center(int atom) {

        atom = atom - 1;
        double xx = x[atom];
        double yy = y[atom];
        double zz = z[atom];
        for (int i = 0; i < N; i++) {
            x[i] = x[i] - xx;
            y[i] = y[i] - yy;
            z[i] = z[i] - zz;
        }
        return this;
    }

    public XYZ times(double factor) {

        for (int i = 0; i < N; i++) {
            x[i] = x[i] * factor;
            y[i] = y[i] * factor;
            z[i] = z[i] * factor;
        }
        return this;
    }

    public XYZ timesX(double factor) {

        for (int i = 0; i < N; i++) {
            x[i] = x[i] * factor;
        }
        return this;
    }

    public XYZ timesY(double factor) {

        for (int i = 0; i < N; i++) {
            y[i] = y[i] * factor;
        }
        return this;
    }

    public XYZ timesZ(double factor) {

        for (int i = 0; i < N; i++) {
            z[i] = z[i] * factor;
        }
        return this;
    }

    public int getCenterAtom() {

        double xSize = Math.abs(Utilities.max(x) - Utilities.min(x));
        double ySize = Math.abs(Utilities.max(y) - Utilities.min(y));
        double zSize = Math.abs(Utilities.max(z) - Utilities.min(z));
        int atom = 0;
        int i = 0;
        double minDistance = Math.sqrt((xSize / 2 - x[i]) * (xSize / 2 - x[i]) + (ySize / 2 - y[i]) * (ySize / 2 - y[i]) + (zSize / 2 - z[i]) * (zSize / 2 - z[i]));;
        for (i = 1; i < N; i++) {

            double distance = Math.sqrt((xSize / 2 - x[i]) * (xSize / 2 - x[i]) + (ySize / 2 - y[i]) * (ySize / 2 - y[i]) + (zSize / 2 - z[i]) * (zSize / 2 - z[i]));
            if (distance < minDistance) {
                minDistance = distance;
                atom = i;
            }

        }
        return atom + 1;
    }

    public XYZ centerX(int atom) {

        double x_atom = x[atom - 1];
        for (int i = 0; i < N; i++) {
            x[i] = x[i] - x_atom;
        }
        return this;
    }

    public XYZ centerX(double xPos) {

        for (int i = 0; i < N; i++) {
            x[i] = x[i] - xPos;
        }
        return this;
    }

    public XYZ centerY(double yPos) {

        for (int i = 0; i < N; i++) {
            y[i] = y[i] - yPos;
        }
        return this;
    }

    public XYZ centerY(int atom) {
        double y_atom = y[atom - 1];
        for (int i = 0; i < N; i++) {
            y[i] = y[i] - y_atom;

        }
        return this;
    }

    public XYZ centerZ(int atom) {
        double z_atom = z[atom - 1];
        for (int i = 0; i < N; i++) {
            z[i] = z[i] - z_atom;

        }
        return this;
    }

    public XYZ centerZ(double zPos) {

        for (int i = 0; i < N; i++) {
            z[i] = z[i] - zPos;
        }
        return this;
    }

    public XYZ rotationalAlignmentOnYAxis(int a, int b) {
        int atomA = a - 1;
        int atomB = b - 1;
        double r12X = x[atomA] - x[atomB];
        double r12Y = y[atomA] - y[atomB];
        double r12Z = z[atomA] - z[atomB];
        double thetaXY = Math.acos(r12X / Math.sqrt(r12X * r12X + r12Y * r12Y));
        double thetaXZ = Math.atan(r12Z / Math.sqrt(r12X * r12X + r12Y * r12Y));

        //Perform rotation 1: about z axis
        for (int i = 0; i < N; i++) {
            double oldx = x[i];
            double oldy = y[i];
            x[i] = oldx * Math.cos(thetaXY) - oldy * Math.sin(thetaXY);
            y[i] = oldx * Math.sin(thetaXY) + oldy * Math.cos(thetaXY);
        }

        //Perform rotation 2: about y axis
        for (int i = 0; i < N; i++) {
            double oldx = x[i];
            double oldz = z[i];
            x[i] = oldx * Math.cos(thetaXZ) + oldz * Math.sin(thetaXZ);
            z[i] = -oldx * Math.sin(thetaXZ) + oldz * Math.cos(thetaXZ);

        }

        double thetaYZ = -Math.atan(z[atomA] / y[atomA]);

        //Perform rotation 3: about x axis
        for (int i = 0; i < N; i++) {
            double oldy = y[i];
            double oldz = z[i];
            y[i] = oldy * Math.cos(thetaYZ) - oldz * Math.sin(thetaYZ);
            z[i] = oldy * Math.sin(thetaYZ) + oldz * Math.cos(thetaYZ);
        }
        return this;
    }

    public XYZ rotationalAlignmentOnZAxis(int a, int b) {
        int atomA = a - 1;
        int atomB = b - 1;
        double r12X = Math.abs(x[atomA] - x[atomB]);
        double r12Y = Math.abs(y[atomA] - y[atomB]);
        double r12Z = Math.abs(z[atomA] - z[atomB]);
        double thetaXY = Math.asin(-r12Y / Math.sqrt(r12X * r12X + r12Y * r12Y));
        double thetaXZ = Math.PI / 2 - Math.atan(-r12Z / Math.sqrt(r12X * r12X + r12Y * r12Y));

        //Perform rotation 1: about z axis
        for (int i = 0; i < N; i++) {
            double oldx = x[i];
            double oldy = y[i];
            x[i] = oldx * Math.cos(thetaXY) - oldy * Math.sin(thetaXY);
            y[i] = oldx * Math.sin(thetaXY) + oldy * Math.cos(thetaXY);
        }

        //Perform rotation 2: about y axis
        for (int i = 0; i < N; i++) {
            double oldx = x[i];
            double oldz = z[i];
            x[i] = oldx * Math.cos(thetaXZ) + oldz * Math.sin(thetaXZ);
            z[i] = -oldx * Math.sin(thetaXZ) + oldz * Math.cos(thetaXZ);

        }
        return this;
    }

    public XYZ rotationalAlignmentAlongZAxis(int a, int b) {
        int atomA = a - 1;
        int atomB = b - 1;
        double r12X = (x[atomA] - x[atomB]);
        double r12Y = (y[atomA] - y[atomB]);
        double r12Z = (z[atomA] - z[atomB]);
        double thetaXY = 0, thetaXZ = 0;
        thetaXY = Math.acos(r12X / Math.sqrt(r12X * r12X + r12Y * r12Y));
        thetaXZ = -Math.PI / 2 + Math.atan(r12Z / Math.sqrt(r12X * r12X + r12Y * r12Y));

        //Perform rotation 1: about z axis
        for (int i = 0; i < N; i++) {
            double oldx = x[i];
            double oldy = y[i];
            x[i] = oldx * Math.cos(thetaXY) - oldy * Math.sin(thetaXY);
            y[i] = oldx * Math.sin(thetaXY) + oldy * Math.cos(thetaXY);
        }

        //Perform rotation 2: about y axis
        for (int i = 0; i < N; i++) {
            double oldx = x[i];
            double oldz = z[i];
            x[i] = oldx * Math.cos(thetaXZ) + oldz * Math.sin(thetaXZ);
            z[i] = -oldx * Math.sin(thetaXZ) + oldz * Math.cos(thetaXZ);

        }
        return this;
    }

    public XYZ remove(int a) {
        a = a - 1;
        String[] tatom = new String[N - 1];
        double[] tx = new double[N - 1];
        double[] ty = new double[N - 1];
        double[] tz = new double[N - 1];
        int j = 0;
        for (int i = 0; i < N; i++) {
            if (i != a) {
                tatom[j] = new String(atom[i]);
                tx[j] = x[i];
                ty[j] = y[i];
                tz[j] = z[i];
                j++;
            }

        }
        x = tx.clone();
        y = ty.clone();
        z = tz.clone();
        atom = tatom.clone();
        N--;
        setup();
        return this;
    }

    public XYZ remove(ArrayList a) {
        ArrayList<String> tatom = new ArrayList();
        ArrayList<Double> tx = new ArrayList();
        ArrayList<Double> ty = new ArrayList();
        ArrayList<Double> tz = new ArrayList();

        int j = 0;
        for (int i = 0; i < N; i++) {
            if (!a.contains(i + 1)) {
                tatom.add(new String(atom[i]));
                tx.add(x[i]);
                ty.add(y[i]);
                tz.add(z[i]);
            }

        }

        atom = new String[tx.size()];

        x = new double[tx.size()];
        y = new double[tx.size()];
        z = new double[tx.size()];

        for (int i = 0; i < tx.size(); i++) {
            x[i] = tx.get(i);
            y[i] = ty.get(i);
            z[i] = tz.get(i);
            atom[i] = tatom.get(i);
        }
        N = atom.length;
        setup();
        return this;
    }

    public XYZ remove(int from, int to) {
        from = from - 1;
        to = to - 1;
        String[] tatom = new String[(to - from + 1)];
        double[] tx = new double[(to - from + 1)];
        double[] ty = new double[(to - from + 1)];
        double[] tz = new double[(to - from + 1)];
        int j = 0;
        for (int i = 0; i < N; i++) {
            if (i <= to && i >= from) {
                tatom[j] = new String(atom[i]);
                tx[j] = x[i];
                ty[j] = y[i];
                tz[j] = z[i];
                j++;
            }
        }
        x = tx.clone();
        y = ty.clone();
        z = tz.clone();
        atom = tatom.clone();
        N = (to - from + 1);
        setup();
        return this;
    }

    public XYZ add(XYZ xyz) {
        String[] tatom = new String[N + xyz.N];
        double[] tx = new double[N + xyz.N];
        double[] ty = new double[N + xyz.N];
        double[] tz = new double[N + xyz.N];

        for (int i = 0; i < N; i++) {
            tatom[i] = new String(atom[i]);
            tx[i] = x[i];
            ty[i] = y[i];
            tz[i] = z[i];
        }

        for (int i = 0; i < xyz.N; i++) {
            tatom[N + i] = new String(xyz.atom[i]);
            tx[N + i] = xyz.x[i];
            ty[N + i] = xyz.y[i];
            tz[N + i] = xyz.z[i];
        }
        N = N + xyz.N;
        atom = new String[N];
        x = new double[N];
        y = new double[N];
        z = new double[N];

        for (int i = 0; i < N; i++) {
            atom[i] = new String(tatom[i]);
            x[i] = tx[i];
            y[i] = ty[i];
            z[i] = tz[i];
        }
        setup();
        return this;
    }

    public XYZ add(String _atom, double _x, double _y, double _z) {

        String[] tatom = new String[N + 1];
        double[] tx = new double[N + 1];
        double[] ty = new double[N + 1];
        double[] tz = new double[N + 1];

        for (int i = 0; i < N; i++) {
            tatom[i] = new String(atom[i]);
            tx[i] = x[i];
            ty[i] = y[i];
            tz[i] = z[i];
        }

        tatom[N] = new String(_atom);
        tx[N] = _x;
        ty[N] = _y;
        tz[N] = _z;
        N = N + 1;
        atom = new String[N];
        x = new double[N];
        y = new double[N];
        z = new double[N];

        for (int i = 0; i < N; i++) {
            atom[i] = new String(tatom[i]);
            x[i] = tx[i];
            y[i] = ty[i];
            z[i] = tz[i];
        }
        setup();
        return this;

    }

    public XYZ add(XYZ a, double _x, double _y, double _z) {

        for (int i = 0; i < a.N; i++) {

        }
        return this;

    }

    public XYZ translate(double tx, double ty, double tz) {
        for (int i = 0; i < N; i++) {
            x[i] += tx;
            y[i] += ty;
            z[i] += tz;
        }
        return this;
    }

    public XYZ makeImages(double[] latticeConstants, int copies, int direction) {
        //center();
        XYZ nxyz = new XYZ(N * copies);
        if (direction == 0)//x axis
        {
            for (int i = 0; i < copies; i++) {
                for (int j = 0; j < N; j++) {
                    nxyz.x[j + i * N] = x[j] + i * latticeConstants[0];
                    nxyz.y[j + i * N] = y[j];
                    nxyz.z[j + i * N] = z[j];
                    nxyz.atom[j + i * N] = new String(atom[j]);
                }
            }
        } else if (direction == 1)//y axis
        {
            for (int i = 0; i < copies; i++) {
                for (int j = 0; j < N; j++) {
                    nxyz.y[j + i * N] = y[j] + i * latticeConstants[1];
                    nxyz.x[j + i * N] = x[j];
                    nxyz.z[j + i * N] = z[j];
                    nxyz.atom[j + i * N] = new String(atom[j]);
                }
            }
        } else if (direction == 2)//z axis
        {
            for (int i = 0; i < copies; i++) {
                for (int j = 0; j < N; j++) {
                    nxyz.z[j + i * N] = z[j] + i * latticeConstants[2];
                    nxyz.y[j + i * N] = y[j];
                    nxyz.x[j + i * N] = x[j];
                    nxyz.atom[j + i * N] = new String(atom[j]);
                }
            }
        }
        setup();
        return nxyz;
    }

    public XYZ makeImages_Hexagonal2D(double latticeConstant, int copies, int direction) {
        XYZ nxyz = new XYZ(N * copies);
        if (direction == 0)//x axis
        {
            for (int i = 0; i < copies; i++) {
                for (int j = 0; j < N; j++) {
                    nxyz.x[j + i * N] = x[j] + i * latticeConstant;
                    nxyz.y[j + i * N] = y[j];
                    nxyz.z[j + i * N] = z[j];
                    nxyz.atom[j + i * N] = new String(atom[j]);
                }
            }
        } else if (direction == 1)//y axis
        {
            for (int i = 0; i < copies; i++) {
                for (int j = 0; j < N; j++) {
                    nxyz.y[j + i * N] = y[j] + i * latticeConstant * Math.cos(Math.PI / 6);
                    nxyz.x[j + i * N] = x[j] + i * latticeConstant * Math.sin(Math.PI / 6);;
                    nxyz.z[j + i * N] = z[j];
                    nxyz.atom[j + i * N] = new String(atom[j]);
                }
            }
        }
        setup();
        return nxyz;
    }

    public XYZ switchAxis(int xNew, int yNew, int zNew) {
        for (int i = 0; i < N; i++) {
            if (xNew == 1 && yNew == 3)//swap y and z
            {
                double temp = y[i];
                y[i] = z[i];
                z[i] = temp;
            } else if (xNew == 2 && yNew == 1)//swap x and y
            {
                double temp = y[i];
                y[i] = x[i];
                x[i] = temp;
            } else if (xNew == 3 && yNew == 2)//swap x and z
            {
                double temp = z[i];
                z[i] = x[i];
                x[i] = temp;
            } else if (xNew == 3 && yNew == 1) {
                double tempx = x[i];
                double tempz = z[i];
                double tempy = y[i];
                z[i] = tempy;
                x[i] = tempz;
                y[i] = tempx;
            }
        }
        return this;
    }

    public XYZ zeroStart(int axis) {
        double minX = Utilities.min(x);
        double minY = Utilities.min(y);
        double minZ = Utilities.min(z);

        for (int i = 0; i < N; i++) {
            if (axis == 1) {
                x[i] = x[i] - minX;
            } else if (axis == 2) {
                y[i] = y[i] - minY;
            } else if (axis == 3) {
                z[i] = z[i] - minZ;
            }
        }
        return this;
    }

    public XYZ rotateZ(double thetaXY) {
        for (int i = 0; i < N; i++) {
            double oldx = x[i];
            double oldy = y[i];
            x[i] = oldx * Math.cos(thetaXY) - oldy * Math.sin(thetaXY);
            y[i] = oldx * Math.sin(thetaXY) + oldy * Math.cos(thetaXY);

        }
        return this;
    }

    public XYZ rotateX(double thetaXY) {

        for (int i = 0; i < N; i++) {
            double oldy = y[i];
            double oldz = z[i];
            y[i] = oldy * Math.cos(thetaXY) - oldz * Math.sin(thetaXY);
            z[i] = oldy * Math.sin(thetaXY) + oldz * Math.cos(thetaXY);
        }
        return this;
    }

    public XYZ rotateY(double thetaXZ) {
        for (int i = 0; i < N; i++) {
            double oldx = x[i];
            double oldz = z[i];
            x[i] = oldx * Math.cos(thetaXZ) + oldz * Math.sin(thetaXZ);
            z[i] = -oldx * Math.sin(thetaXZ) + oldz * Math.cos(thetaXZ);

        }
        return this;
    }

    public static XYZ constructElectrodeScattering(XYZ left, XYZ scattering, XYZ right, double spacing) {
        XYZ result = new XYZ(scattering.N + left.N + right.N);

        left.zeroStart(3);
        right.zeroStart(3);
        scattering.zeroStart(3);

        for (int i = 0; i < left.N; i++) {
            result.x[i] = left.x[i];
            result.y[i] = left.y[i];
            result.z[i] = left.z[i];
            result.atom[i] = new String(left.atom[i]);
        }

        double LeftZ = Utilities.max(left.z) - Utilities.min(left.z);

        for (int i = 0; i < scattering.N; i++) {
            result.x[i + left.N] = scattering.x[i];
            result.y[i + left.N] = scattering.y[i];
            result.z[i + left.N] = LeftZ + spacing + scattering.z[i];
            result.atom[i + left.N] = new String(scattering.atom[i]);
        }

        double maxLeftAndScatteringZ = LeftZ + spacing + Utilities.max(scattering.z) - Utilities.min(scattering.z);

        for (int i = 0; i < right.N; i++) {
            result.x[i + left.N + scattering.N] = right.x[i];
            result.y[i + left.N + scattering.N] = right.y[i];
            result.z[i + left.N + scattering.N] = maxLeftAndScatteringZ + spacing + right.z[i];
            result.atom[i + left.N + scattering.N] = new String(right.atom[i]);
        }
        return result;
    }

    public static XYZ constructElectrodeScattering(XYZ left, XYZ scattering, XYZ right, int scattL, int toLeft, int scattR, int toRight) {
        //Removing the two atoms that were placed as placeholders to be substituted with those 
        //terminal atoms on the two electrodes.
        XYZ result = new XYZ(scattering.N - 2 + left.N + right.N);

        left.zeroStart(3);
        right.zeroStart(3);
        scattering.zeroStart(3);

        left.centerX(toLeft);
        left.centerY(toLeft);
        right.centerX(toRight);
        right.centerY(toRight);
        scattering.centerX(scattL);
        scattering.centerY(scattL);

        for (int i = 0; i < left.N; i++) {
            result.x[i] = left.x[i];
            result.y[i] = left.y[i];
            result.z[i] = left.z[i];
            result.atom[i] = new String(left.atom[i]);
        }

        double LeftZ = Utilities.max(left.z) - Utilities.min(left.z);

        int resultCounter = 0;
        for (int i = 0; i < scattering.N; i++) {
            if (i + 1 == scattL || i + 1 == scattR) {
                System.out.println(scattering.atom[i]);
                continue;
            }
            result.x[resultCounter + left.N] = scattering.x[i];
            result.y[resultCounter + left.N] = scattering.y[i];
            result.z[resultCounter + left.N] = LeftZ + scattering.z[i];
            result.atom[resultCounter + left.N] = new String(scattering.atom[i]);
            resultCounter++;

        }

        double maxLeftAndScatteringZ = LeftZ + Utilities.max(scattering.z) - Utilities.min(scattering.z);

        for (int i = 0; i < right.N; i++) {
            result.x[i + left.N + scattering.N - 2] = right.x[i];
            result.y[i + left.N + scattering.N - 2] = right.y[i];
            result.z[i + left.N + scattering.N - 2] = maxLeftAndScatteringZ + right.z[i];
            result.atom[i + left.N + scattering.N - 2] = new String(right.atom[i]);
        }
        return result;
    }

    public XYZ translatePeriodic(double[] latticeConstant, double[] translation) {
        double minX = Utilities.min(x);
        double minY = Utilities.min(y);
        double minZ = Utilities.min(z);

        for (int i = 0; i < N; i++) {
            double trans = x[i] - translation[0];
            if (trans < minX) {
                x[i] = x[i] - translation[0] + latticeConstant[0];
            } else {
                x[i] = x[i] - translation[0];
            }
            trans = y[i] - translation[1];
            if (trans < minY) {
                y[i] = y[i] - translation[1] + latticeConstant[1];
            } else {
                y[i] = y[i] - translation[1];
            }
            trans = z[i] - translation[2];
            if (trans < minZ) {
                z[i] = z[i] - translation[2] + latticeConstant[2];
            } else {
                z[i] = z[i] - translation[2];
            }
        }
        return this;
    }

    public Matrix getCoulombMatrix(String[] atomsArray) {
        int n = atomsArray.length;
        Matrix c = new Matrix(n, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    c.set(i, i, Math.pow(AtomicData.getMassNumber(atom[i]), 2.4));
                } else {
                    c.set(i, j, AtomicData.getMassNumber(atom[i]) * AtomicData.getMassNumber(atom[j]) / distance(i + 1, j + 1));
                }
            }
        }

        for (int i = N; i < n; i++) {
            for (int j = N; j < n; j++) {
                c.set(i, j, 0);
            }
        }

        return c;
    }

    public int getNumAtoms(String label) {
        int t = 0;
        for (int i = 0; i < N; i++) {
            if (label.equals(atom[i])) {
                t++;
            }
        }
        return t;
    }

    public int getNumValenceElectrons() {
        int valence = 0;
        for (int i = 0; i < N; i++) {
            valence += AtomicData.getValence(atom[i]);
        }
        return valence;
    }

    public void sortByLabel() {
        for (int i = 0; i < N - 1; i++) {
            for (int j = i; j < N; j++) {
                if (atom[i].compareTo(atom[j]) < 0) {
                    String t = atom[i];
                    atom[i] = atom[j];
                    atom[j] = t;

                    double tx = x[i];
                    x[i] = x[j];
                    x[j] = tx;

                    tx = y[i];
                    y[i] = y[j];
                    y[j] = tx;

                    tx = z[i];
                    z[i] = z[j];
                    z[j] = tx;
                }
            }
        }
    }

    public double[] getContainingBox(double distanceFromBoundary) {
        double xSize = Math.abs(Utilities.max(x) - Utilities.min(x));
        double ySize = Math.abs(Utilities.max(y) - Utilities.min(y));
        double zSize = Math.abs(Utilities.max(z) - Utilities.min(z));

        return new double[]{xSize + distanceFromBoundary * 2, ySize + distanceFromBoundary * 2, zSize + distanceFromBoundary * 2};
    }

    public void addRandomly(String atom, double latticeConstant) {
        double xx = 0;
        double yy = 0;
        double zz = 0;

        boolean done = false;
        while (!done) {
            boolean failed = false;
            for (int i = 0; i < N; i++) {
                if (distance(i + 1, new double[]{xx, yy, zz}) < 1.5) {
                    failed = true;
                    break;
                }
            }
            if (!failed) {
                done = true;
                add(atom, xx, yy, zz);
                return;
            } else {
                xx = latticeConstant * Math.random();
                yy = latticeConstant * Math.random();
                zz = latticeConstant * Math.random();
            }
        }
    }

    public void addRandomly(String atom, int num, double from, double to) {
        double xx = 0;
        double yy = 0;
        double zz = 0;
        center();
        double xlatticeConstant = Math.abs(Utilities.max(x) - Utilities.min(x)) + to;
        double ylatticeConstant = Math.abs(Utilities.max(y) - Utilities.min(y)) + to;
        double zlatticeConstant = Math.abs(Utilities.max(z) - Utilities.min(z)) + to;
        XYZ added = new XYZ(0);
        for (int atoms = 0; atoms < num; atoms++) {
            boolean done = false;

            while (!done) {
                xx = xlatticeConstant * Math.random();
                yy = ylatticeConstant * Math.random();
                zz = zlatticeConstant * Math.random();

                boolean failed = false;
                for (int i = 0; i < N; i++) {
                    if (distance(i + 1, new double[]{xx, yy, zz}) < from || distance(i + 1, new double[]{xx, yy, zz}) > to) {
                        failed = true;
                        break;
                    }
                }
                if (!failed) {
                    System.out.println("Added at " + xx + "," + yy + "," + zz);
                    done = true;
                    added.add(atom, xx, yy, zz);
                    return;
                }
            }

        }
        add(added);
    }

    public void addRandomlyBruteForce(String atom, ArrayList addedAtomX, ArrayList addedAtomY, ArrayList addedAtomZ, int num, double from, double to) {

        center();
        double xlatticeConstant = Math.abs(Utilities.max(x) - Utilities.min(x)) + to + 5;
        double ylatticeConstant = Math.abs(Utilities.max(y) - Utilities.min(y)) + to + 5;
        double zlatticeConstant = Math.abs(Utilities.max(z) - Utilities.min(z)) + to + 5;
        double xStep = xlatticeConstant / 100;
        double yStep = ylatticeConstant / 100;
        double zStep = zlatticeConstant / 100;
        boolean[][][] box = new boolean[100][100][100];
        ArrayList arrayOfXPoints = new ArrayList();
        ArrayList arrayOfYPoints = new ArrayList();
        ArrayList arrayOfZPoints = new ArrayList();
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 100; j++) {
                for (int k = 0; k < 100; k++) {
                    box[i][j][k] = false;
                }
            }
        }
        for (int a = 0; a < N; a++) {
            //label all points that live on the union of spheres surrounding atoms
            for (int i = 0; i < 100; i++) {
                for (int j = 0; j < 100; j++) {
                    for (int k = 0; k < 100; k++) {
                        double xPos = (i - 50) * xStep;
                        double yPos = (j - 50) * yStep;
                        double zPos = (k - 50) * zStep;

                        if ((distance(a + 1, new double[]{xPos, yPos, zPos}) > from && distance(a + 1, new double[]{xPos, yPos, zPos}) < to)) {
                            box[i][j][k] = true;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 100; j++) {
                for (int k = 0; k < 100; k++) {
                    if (box[i][j][k]) {
                        boolean ok = true;
                        double xPos = (i - 50) * xStep;
                        double yPos = (j - 50) * yStep;
                        double zPos = (k - 50) * zStep;
                        for (int a = 0; a < N; a++) {
                            if (distance(a + 1, new double[]{xPos, yPos, zPos}) < from) {
                                ok = false;
                                break;
                            }
                        }
                        if (ok) {
                            arrayOfXPoints.add(xPos);
                            arrayOfYPoints.add(yPos);
                            arrayOfZPoints.add(zPos);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < num; i++) {
            int index = (int) (Math.random() * arrayOfXPoints.size());
            double xPos = (Double) arrayOfXPoints.get((index));
            double yPos = (Double) arrayOfYPoints.get((index));
            double zPos = (Double) arrayOfZPoints.get((index));
            boolean bad = false;
            for (int k = 0; k < addedAtomX.size(); k++) {
                double addedX = (double) addedAtomX.get(k);
                double addedY = (double) addedAtomY.get(k);
                double addedZ = (double) addedAtomZ.get(k);

                if (Math.sqrt(Math.pow(addedX - xPos, 2) + Math.pow(addedY - yPos, 2) + Math.pow(addedZ - zPos, 2)) < 3) {
                    bad = true;
                    break;
                }
            }
            if (!bad) {
                add(atom, xPos, yPos, zPos);
                addedAtomX.add(xPos);
                addedAtomY.add(yPos);
                addedAtomZ.add(zPos);
            } else {
                i = i - 1;
            }
        }
    }

    public void addRandomly(String atom, double bond, double a, double b, double c) {
        double xx = 0;
        double yy = 0;
        double zz = 0;

        boolean done = false;
        while (!done) {
            boolean failed = false;
            for (int i = 0; i < N; i++) {
                if (distance(i + 1, new double[]{xx, yy, zz}) < bond) {
                    failed = true;
                    break;
                }
            }
            if (!failed) {
                done = true;
                add(atom, xx, yy, zz);
                return;
            } else {
                xx = a * Math.random();
                yy = b * Math.random();
                zz = c * Math.random();
            }
        }
    }

    public boolean addRandomlyUntilFull(String atom, double a, double b, double c) {
        double xx = 0;
        double yy = 0;
        double zz = 0;

        boolean done = false;
        for (int t = 0; t < 10000; t++) {
            boolean failed = false;
            for (int i = 0; i < N; i++) {
                if (distance(i + 1, new double[]{xx, yy, zz}) < 1.5) {
                    failed = true;
                    break;
                }
            }
            if (!failed) {
                done = true;
                add(atom, xx, yy, zz);
                return true;
            } else {
                xx = a * Math.random();
                yy = b * Math.random();
                zz = c * Math.random();
            }
        }
        return done;
    }

    public void addRandomly(XYZ molecule, double bond, double a, double b, double c) {
        double xx = 0;
        double yy = 0;
        double zz = 0;

        boolean done = false;
        while (!done) {
            boolean failed = false;
            for (int i = 0; i < N; i++) {
                if (distance(i + 1, molecule) < bond) {
                    failed = true;
                    break;
                }
            }
            if (!failed) {
                done = true;
                add(molecule);
                return;
            } else {
                xx = a * Math.random();
                yy = b * Math.random();
                zz = c * Math.random();

                double xAngle = Math.PI * 2 * Math.random();
                double yAngle = Math.PI * 2 * Math.random();
                double zAngle = Math.PI * 2 * Math.random();
                molecule.center();
                molecule.rotateX(xAngle);
                molecule.rotateY(yAngle);
                molecule.rotateZ(zAngle);

                molecule.translate(xx, yy, zz);
            }
        }
    }

    public void cut(double centerXp, double centerYp, double radius) {
        centerX(centerXp);
        centerY(centerYp);
        ArrayList a = new ArrayList();
        for (int i = 1; i <= x.length; i++) {
            if (distance2D(i, new double[]{0, 0}) > radius) {
                a.add(i);
            }
        }
        remove(a);
    }

    /*
    When you want to sort your current xyz to resemble that of an inputXYZ
    For example: you want to use the coordinates file from VASP with phonon calculations in SIESTA
     */
    public XYZ sortAccordingTo(XYZ inputXYZ, int baseAtomHere, int baseAtomThere) {
        ArrayList<String> tatom = new ArrayList();
        ArrayList<Double> tx = new ArrayList();
        ArrayList<Double> ty = new ArrayList();
        ArrayList<Double> tz = new ArrayList();

        inputXYZ.centerX(baseAtomThere);
        inputXYZ.centerY(baseAtomThere);

        centerX(baseAtomHere);
        centerY(baseAtomHere);

        System.out.println(N + ", " + inputXYZ.N);

        for (int j = 0; j < inputXYZ.N; j++) {
            for (int i = 0; i < N; i++) {
                if (distance(i + 1, new double[]{inputXYZ.x[j], inputXYZ.y[j], inputXYZ.z[j]}) < 1) {
                    tatom.add(atom[i]);
                    tx.add(x[i]);
                    ty.add(y[i]);
                    tz.add(z[i]);
                    break;
                }
            }
        }
        atom = new String[tx.size()];

        x = new double[tx.size()];
        y = new double[tx.size()];
        z = new double[tx.size()];

        for (int i = 0; i < tx.size(); i++) {
            x[i] = tx.get(i);
            y[i] = ty.get(i);
            z[i] = tz.get(i);
            atom[i] = tatom.get(i);
        }
        N = atom.length;
        return this;
    }

    public XYZ createRipple(int axis, double length, double amplitude) {
        if (axis == 2)//y-axis
        {
            for (int i = 0; i < N; i++) {
                z[i] += z[i] + amplitude * Math.sin(y[i] * 2 * Math.PI / length);
            }
        } else if (axis == 1)//y-axis
        {
            for (int i = 0; i < N; i++) {
                z[i] += z[i] + amplitude * Math.sin(x[i] * 2 * Math.PI / length);
            }
        }

        return this;
    }

    /*
    Processes the file, which is extracted from OUTCAR using:
    grep -A165 "POSITION                                       TOTAL-FORCE (eV/Angst)" OUTCAR-1 > p_1
    
    Here, -A165 is for a system with 164 atoms
     */
    public static ArrayList<XYZ> getXYZFromVASPMD(String path, String file, String[] poscarLabels, int[] poscarCounts) {
        ArrayList<XYZ> xyzList = new ArrayList<>();
        BufferedReader in = null;
        try {
            ArrayList<String> atoms = new ArrayList();
            in = new BufferedReader(new FileReader(new File(path + file)));
            int numAtoms = 0;
            for (int i = 0; i < poscarCounts.length; i++) {
                numAtoms += poscarCounts[i];
                for (int j = 0; j < poscarCounts[i]; j++) {
                    atoms.add(poscarLabels[i]);
                }
            }

            int iter = 0;
            String line = "";
            line = in.readLine();
            while (line != null) {

                XYZ xyz = new XYZ(numAtoms);
                for (int i = 0; i < numAtoms; i++) {

                    xyz.atom[i] = atoms.get(i);

                    String[] parts = line.split("\\s+");
                    xyz.x[i] = Double.parseDouble(parts[1]);
                    xyz.y[i] = Double.parseDouble(parts[2]);
                    xyz.z[i] = Double.parseDouble(parts[3]);
                    line = in.readLine();
                }

                xyzList.add(xyz);
                iter++;
            }
        } catch (IOException e) {
            e.printStackTrace();

        } finally {
            if (in != null) {
                try {
                    in.close();
                } catch (IOException ioe) {
                    ioe.printStackTrace();
                }
            }
        }
        return xyzList;
    }
}
