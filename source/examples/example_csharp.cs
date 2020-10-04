using System;
using System.Runtime.InteropServices;
using Xatlas;

namespace XatlasExample
{
    class Program
    {
        static void Main(string[] args)
        {
            xatlas.ProgressCategory pc = xatlas.ProgressCategory.ComputeCharts;
            IntPtr namePtr = xatlas.ProgressCategoryString(pc);
            string name = Marshal.PtrToStringAnsi(namePtr);
            Console.WriteLine(name);
            Console.ReadKey();
        }
    }
}
