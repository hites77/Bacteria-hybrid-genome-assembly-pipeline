nextflow.enable.dsl=2

/**
 * Get longest common directory of a list of files.
 */
def getDirectory(fileList) {
    // make paths absolute
    for (int i=0; i < fileList.size(); i++) {
        fileList[i] = fileList[i].toAbsolutePath()
    }

    // try to find longest common directory
    def directory = fileList[0].isDirectory() ? fileList[0] : file(fileList[0].parent)
    boolean continueFlag = true
    while (continueFlag) {
        continueFlag = false
        for (int i=0; i < fileList.size(); i++) {
            if (fileList[i] != directory) {
                continueFlag = true
                if (fileList[i].toString().length() >= directory.toString().length()) {
                    fileList[i] = file(fileList[i].parent)
                }

                if (fileList[i].toString().length() < directory.toString().length()) {
                    directory = fileList[i]
                }
            }
        }
    }

    return directory
}

/**
 * Transforms a channel of lists of files to a channel of directories by applying getDirectory to each list.
 */
def mapToDirectory(fileListChan) {
    return fileListChan.map { getDirectory(it) }
}

def getCondaEnv(envName) {
    return params.condaEnvsDir + '/' + envName
}
