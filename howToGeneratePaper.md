1. Download docker container:  
`> docker pull openjournals/paperdraft`

2. Run container and share the current directory with it (this repository) (${PWD} is powershell syntax for present working directory)
`> docker run --rm -v ${PWD}:/data --env JOURNAL=joss openjournals/paperdraft`