
function rename_files(folder::String, old_substring::String, new_substring::String)
    files = readdir(folder)

    for file in files
        old_path = joinpath(folder, file)
        
        if isfile(old_path)
            new_name = replace(file, old_substring => new_substring)
            new_path = joinpath(folder, new_name)

            if new_name != file
                println("Renaming '$file' to '$new_name'")
                mv(old_path, new_path)
            end
        end
    end
end

folder = "./solutions/problem_P13/DRS_FP_r0"
old_substring = "idx"
new_substring = "d_100_idx"

rename_files(folder, old_substring, new_substring)