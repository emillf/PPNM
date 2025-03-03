public class Example {
  public static void Main(string[] args) {
    try {
      // Some code that might throw an exception
		int zero=0;
      int result = 10 / zero; // Example: Division by zero
    } catch (System.DivideByZeroException ex) {
      System.Console.Error.WriteLine("An error occurred: " + ex.Message);
      System.Environment.ExitCode = 1; // Indicate an error to the operating system
    }
  }
}
