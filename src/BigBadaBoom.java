/*
 * Algorithms Assignment 2, Problem 3
 * Andrew Stagg
 * Adam Chapman
 * Brad LaVigne
 **/
import java.io.*;
import java.util.*;

public class BigBadaBoom {
   public static void main(String[] args) throws Exception {
      Compound[] cees;
      cees = read_data("input3.txt");
      // cees = create_data(10000000, 10, 32768);

      Long t = System.currentTimeMillis();
      Arrays.sort(cees);
      t = System.currentTimeMillis() - t;
      // System.out.println("Sorted in " + t + " ms");

      t = System.currentTimeMillis();
      boom_rec(cees, 0, cees.length - 1);
      t = (System.currentTimeMillis() - t);

      System.out.println(Compound.mu.id);
      System.out.println(Compound.nu.id);
   }

   /* Divide the sorted array at the median and recurse */
   static void boom_rec(Compound[] cees, int left, int right) {
      if (left == right) {
      } else if (right - left == 1) {
         cees[left].react(cees[right]);
      } else {
         int mid = (left + right) / 2;
         boom_rec(cees, left, mid);
         boom_rec(cees, mid + 1, right);

         // for very short subarrays, use cross_C due to lower overhead
         if (right - left < 8) cross_C(cees, left, mid, right);
         else cross_N(cees, left, mid, right);
      }
   }

   /* Search crossover region based on existing Carbon sort */
   static void cross_C(Compound[] cees, int left, int mid, int right) {
      int lo = (int) (cees[mid].C - Compound.E + 0.99);
      int hi = (int) (cees[mid].C + Compound.E);

      // For each compound l on left, react with each compound r on right
      // pre-sorted by C, so early termination on right when r.C-l.C>=E
      for (int l = mid; l >= left && cees[l].C >= lo; l--) {
         hi = (int) (cees[l].C + Compound.E); // max r.C to match with this l.C
         for (int r = mid + 1; r <= right && cees[r].C <= hi; r++) {
            // Short-circuit if N or O out of range. Not much help when Cross_N
            // is in play, but gives a 2x improvement when only using Cross_C
            if (!(cees[r].N - cees[l].N > Compound.E || cees[r].O - cees[l].O > Compound.E))
               cees[l].react(cees[r]);
         }
      }
   }

   /*
    * Search crossover region based on new sort by Nitrogen
    * This provides additional sparsity compared to the existing Carbon sort,
    * thus requiring fewer cross-checks
    * 
    * The original sort has N as second priority and there seems to be a benefit
    * from this that carries over into the secondary sort, since the data is
    * already partially sorted by N. Experimentally, I changed to an O sort here
    * and performance suffered.
    */
   static void cross_N(Compound[] cees, int left, int mid, int right) {
      int lo = (int) (cees[mid].C - Compound.E + 0.99);
      int l = mid;
      while (l > left && cees[l - 1].C >= lo)
         l--;// iterate l to lower bound
      int hi = (int) (cees[mid].C + Compound.E);
      int r = mid;
      while (r < right && cees[r + 1].C <= hi)
         r++;// iterate r to upper bound

      // Create distinct arrays for each side of the dividing plane
      Compound[] front = Arrays.copyOfRange(cees, l, mid);
      Compound[] back = Arrays.copyOfRange(cees, mid == r ? mid : mid + 1, r);

      // Sort each side by N
      Compound.secondary_sort = true;
      Arrays.sort(front);
      Arrays.sort(back);
      Compound.secondary_sort = false;

      // Because the compounds are spread over a wider range in N than in C, the
      // early termination clause becomes so more effective that adding any sort
      // of test against C or O reduces performance by a factor similar to that
      // by which such tests benefit Cross_C.
      for (Compound c : front)
         for (Compound d : back)
            if (d.N - c.N > (int) Compound.E) break;
            else c.react(d);
   }

   /* Read data from file */
   static Compound[] read_data(String file) {
      List<Compound> list = new ArrayList<Compound>(1000000);
      BufferedReader in;
      try {
         in = new BufferedReader(new FileReader(file));
         while (in.ready()) {
            String[] t = in.readLine().split(",");
            list.add(new Compound(t));
         }
         in.close();
      } catch (IOException e) {
         e.printStackTrace();
      }
      return list.toArray(new Compound[0]);
   }

   /*
    * Generate randomized data such that good solutions are uncommon and
    * "perfect" solutions (E=1.0) are unlikely
    */
   static Compound[] create_data(int n, int m, int v) {
      LinkedHashSet<Compound> set = new LinkedHashSet<Compound>(n);
      String[] t = new String[4];

      for (int i = 1; i <= m; i++) {
         while (set.size() < Math.pow(n, 1.0 / m * i)) {
            t[1] = "" + (int) (Math.random() * v / i) * i;
            t[2] = "" + (int) (Math.random() * v / i) * i;
            t[3] = "" + (int) (Math.random() * v / i) * i;
            t[0] = "" + (set.size() + 1);
            set.add(new Compound(t));
         }
      }
      Compound[] cees = set.toArray(new Compound[0]);
      return cees;
   }

   /* Write data to file */
   void write_data(Compound[] cees, String file) {
      try {
         FileWriter f = new FileWriter(file, false);
         for (Compound c : cees)
            f.write(c.id + "," + c.C + "," + c.N + "," + c.O + "\n");
         f.close();
      } catch (IOException e) {
         e.printStackTrace();
      }
   }
}

class Compound implements Comparable<Compound> {
   static boolean secondary_sort = false;
   static long count = 0; // calls to react
   static Compound mu = null; // global best
   static Compound nu = null; // global best
   static double E = Double.MAX_VALUE; // global best
   int id;
   int C;
   int N;
   int O;

   Compound(String[] t) {
      this.id = Integer.parseInt(t[0]);
      this.C = Integer.parseInt(t[1]);
      this.N = Integer.parseInt(t[2]);
      this.O = Integer.parseInt(t[3]);
      mu = nu = this;
   }

   /* Compute reactivity of two compounds */
   void react(Compound c) {
      // count++;
      double e = this.equals(c) ? Double.MAX_VALUE : Math.sqrt((C - c.C)
            * (C - c.C) + (N - c.N) * (N - c.N) + (O - c.O) * (O - c.O));
      if (e < E) {
         mu = this;
         nu = c;
         E = e;
      }
   }

   /* Thinking of my dear friend Dr. Krohn... */
   public int compareTo(Compound c) {
      return secondary_sort ? N < c.N ? -1 : N > c.N ? 1 : C < c.C ? -1
            : C > c.C ? 1 : O < c.O ? -1 : O > c.O ? 1 : 0 : C < c.C ? -1
            : C > c.C ? 1 : N < c.N ? -1 : N > c.N ? 1 : O < c.O ? -1
                  : O > c.O ? 1 : 0;
   }

   public boolean equals(Compound that) {
      return C == that.C && N == that.N && O == that.O;
   }

   public String toString() {
      return "id " + id + "\tC " + C + "\tN " + N + "\tO " + O;
   }
}