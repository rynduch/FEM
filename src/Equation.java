public class Equation {
  double[] t;
  double[][] temp;

  public Equation(Aggregation A) {
    setT(A);
    setTemp(A);
  }

  public void display() {
    System.out.println();
    System.out.println("----------------------------------------------------------------------------------------------------");
    System.out.println("                                           EQUATIONS");
    System.out.println("----------------------------------------------------------------------------------------------------");
    System.out.println();
    displayT();
    displayTemp();
  }

  public void setT(Aggregation A) {
    t = new double[A.HHbc_G.length];
    t = solveEquation(A.HHbc_G, A.P_G);
  }

  public void displayT() {
    System.out.println("[H]{t} + [P] = 0");
    System.out.println("t = {");
    for (double v : t) {
      System.out.print("\t" + v);
      System.out.println();
    }
    System.out.println("}");
  }

  public void displayTemp() {
    System.out.println("([H]+[C]/Δτ){t1}−([C]/Δτ){t0}+P=0");
    for (int i = 0; i < temp.length; i++) {
      double t_max = temp[i][0];
      double t_min = temp[i][0];
      System.out.printf("t%d = {", i + 1);
      System.out.println();

      for (int j = 1; j < t.length; j++) {
        if (temp[i][j] > t_max)
          t_max = temp[i][j];
        else if (temp[i][j] < t_min)
          t_min = temp[i][j];
      }
      System.out.print("\tmin: " + t_min);
      System.out.print("\tmax: " + t_max);
      System.out.println();
      System.out.println("}");
    }
  }

  public void setTemp(Aggregation A) {
    temp = new double[A.aGrid.simulationTime / A.aGrid.simulationStepTime][A.CdTT0P_G.length];
    temp[0] = solveEquation(A.HCdT_G, A.CdTT0P_G);
    double[][] CdTTP = new double[A.CdTT0P_G.length][A.HCdT_G.length];
    for (int w = 0; w < CdTTP.length; w++) {
      for (int k = 0; k < CdTTP[0].length; k++) {
        CdTTP[w][k] = 0;
      }
    }
    for (int t = 1; t < temp.length; t++) {
      for (int w = 0; w < CdTTP[0].length; w++) {
        for (int k = 0; k < CdTTP[0].length; k++) {
          CdTTP[t][w] += A.CdT_G[w][k] * temp[t - 1][k];
        }
        CdTTP[t][w] += A.P_G[w];
      }
      temp[t] = solveEquation(A.HCdT_G, CdTTP[t]);
    }
  }

  public static double[] solveEquation(double[][] A, double[] B) {
    int n = A.length;
    double factor;

    double[][] matrix = new double[n][n + 1];
    for (int i = 0; i < n; i++) {
      System.arraycopy(A[i], 0, matrix[i], 0, n);
      matrix[i][n] = B[i];
    }

    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        factor = matrix[j][i] / matrix[i][i];
        for (int k = i; k <= n; k++) {
          matrix[j][k] -= factor * matrix[i][k];
        }
      }
    }

    double[] result = new double[n];
    for (int i = n - 1; i >= 0; i--) {
      result[i] = matrix[i][n];
      for (int j = i + 1; j < n; j++) {
        result[i] -= matrix[i][j] * result[j];
      }
      result[i] /= matrix[i][i];
    }

    return result;
  }
}
