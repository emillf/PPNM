public class Example {
  public static void Main(string[] args) {
    System.Console.Write("Enter your name: "); // Write without a newline
    string name = System.Console.ReadLine(); // Read a line

    System.Console.Write("Enter your age: ");
    string ageString = System.Console.ReadLine();
    int age = int.Parse(ageString); // Convert string to integer (error handling omitted for brevity)

    System.Console.WriteLine($"Hello, {name}! You are {age} years old.");
  }
}
