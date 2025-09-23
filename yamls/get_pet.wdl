#
version 1.0

workflow WhoAmI {
  call GetPetAccount
}

task GetPetAccount {
  command <<<
    set -euo pipefail
    echo "Active service account:"
    gcloud config get-value account > pet_account.txt
  >>>

  output {
    File pet_account = "pet_account.txt"
  }

  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"   # lightweight docker with curl included
    memory: "1G"
    preemptible: 1
  }
}

