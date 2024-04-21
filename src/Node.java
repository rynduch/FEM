public class Node {
  double x;
  double y;
  public int BC = 0;

  public Node() {
    this.x = 0;
    this.y = 0;
    this.BC = 0;
  }

  public Node(double x, double y) {
    this.x = x;
    this.y = y;
  }

  public void setBC(int bc) {
    this.BC = bc;
  }
}
