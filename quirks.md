# Known code quality and developer environment issues

## Quirks

* I had to specify x64 architecture.
  * The compiler thinks QDK is only compiled for that architecture.
  * Is this true, or bad settings on my part?
* I have no idea what I'm doing with the .csproj files.
* From time to time, the build hangs.
  * I `taskkill /f /im dotnet.exe` to get control back (ctrl-c from command-line build does nothing). Running the exact same build command then produces a successful build, sometimes after a few more taskkills.
  * VS Code's compiler integration magic might be the culprit somehow.
  * Once the build succeeds, multiple command-line builds in repeated succession result in repeated success.
  * taskkill often, but not always, kills multiple `dotnet.exe` instances. Why are there multiple ones running?
  * Will try clean reinstall of qdk.
* My TermsDictionary class seems odd.
  * Why do I need to instantiate a comparer object for things to work? Why not just a (static class) function?
  * Are my ImmutableArray usages correct?
* Everything (tested) is public.
* I have tried to follow QDK documentation conventions but have at times felt the need to invent new headings which may not be parsed correctly by automatic documentation software.

## Things I know I've done badly

* I may still need to refactor code to impose standard capitalization conventions.
* I need to understand standard naming conventions w.r.t. c#'s "everything is in a class" rule. The "FSTools" class is almost certainly misnamed.

## To Do

* Produce an exception or output flag when the swap network can't locally evaluate every term in the Hamiltonian, instead of silently doing partial evaluation.
* Complete comment documentation for recent parts of the code
* Build a demo that evaluates 2D Hubbard
* Refactor library portion of code as a QDK enhancement so that everything can be driven from Q#
