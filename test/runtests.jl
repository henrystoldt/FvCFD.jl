# Runs tests from any file which which starts with "test_" and ends with ".jl"

# First assume the file is running from ./test (this occurs when running Pkg.test)
files = readdir(".")
pathPrefix = ""

# Check if we are actually running from the main directory, and read files from ./test instead
    # this occurs when running include("runtests.jl") from a julia instance in the main directory
if "src" in files
    files = readdir("./test")
    pathPrefix = "./test/"
end

# Look for and run test files
for file in files
    if lastindex(file) > 5 && file[1:5] == "test_" && file[end-2:end] == ".jl"
        println("Running Tests from: test/$file")
        include("$file")
        println()
    end
end