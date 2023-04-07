# Known code quality and developer environment issues

## Quirks

* I had to specify x64 architecture.
  * The compiler thinks QDK is only compiled for that architecture.
  * Is this true, or bad settings on my part?
* I have no idea what I'm doing with the .csproj files.
* I have tried to follow QDK documentation conventions but have at times felt the need to invent new headings which may not be parsed correctly by automatic documentation software.

## Things I know I've done badly

* Is there a naming convention that improves on "FSTools"?

## Desiderata

* Produce an exception or output flag when the swap network can't locally evaluate every term in the Hamiltonian, instead of silently doing partial evaluation.
* Complete comment documentation for recent parts of the code
* Refactor library portion of code as a QDK enhancement so that everything can be driven from Q#
