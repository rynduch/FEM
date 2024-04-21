public class Main {

  public static void main(String[] args) {

    int n = 2;     // Gauss n-point formula, n∈{2,3,4}
    UnivElement ue = new UnivElement(n);
    ue.display();
    Grid g = new Grid();
    g.importData("example.txt");
    g.displayData();
    g.count(ue);
    g.display();
    Aggregation a = new Aggregation(g);
    a.display();
    Equation e = new Equation(a);
    e.display();

  }
}