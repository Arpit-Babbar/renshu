using PkgTemplates
t = Template(;user="arpit-babbar", authors = "Arpit Babbar",
	      julia = VERSION, plugins = [
					  Git,
					  GitHubActions(; x86 = true, windows = true, osx = true),
				          Codecov(),
					  Documenter{GitHubActions}()
				         ]
	     )
t("SpectralMethodsTutorials")

# After this, make a new repository and then set its origin
# $ git remote rm origin
# $ git remote add origin git@github.com:arpit-babbar/repo_name
# $ git config master.remote origin
# $ git config master.merge refs/heads/master
# source - https://gist.github.com/DianaEromosele/fa228f6f6099a8996d3cb891109ab975
