nextflow.enable.dsl=2

/**
 * Returns a closure to be used with publishDir's saveAs parameter which ensures
 * .command.sh, .command.log and .command.sh are be published to params.oudir + params.o_pubdir.
 *
 * @param o_pubdir: params.o.processName eg. process.o.cleanShortReads
 *
 */
def makeNextflowLogClosure(o_pubdir) {
    return { // it = file name
        if (it == '.exitcode' || it == '.command.log' || it == '.command.sh' ) {
            return params.outdir + o_pubdir + 'nextflow' + it
        } else {
            return it
        }
    }
}

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
    boolean continueFlag = false
    while (true) {
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
        if (!continueFlag) {
            break
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
