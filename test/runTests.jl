# Runs tests from any file which which starts with "test_" and ends with ".jl"
files = readdir("./Test")

for file in files
    if lastindex(file) > 5 && file[1:5] == "test_" && file[end-2:end] == ".jl"
        println("Running Tests from: test/$file")
        include("$file")
        println()
    end
end
