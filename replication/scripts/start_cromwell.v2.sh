# Create cromwell config file
if ! [ -e /home/jupyter/.cromwell ]; then
  mkdir /home/jupyter/.cromwell
fi
if [ -e /home/jupyter/.cromwell/cromwell.conf ]; then
  rm /home/jupyter/.cromwell/cromwell.conf
  rm /home/jupyter/.cromwell/cromwell.override.conf
fi

#wb cromwell generate-config \
#  --google-bucket-name=$WORKSPACE_BUCKET \
#  --dir=/home/jupyter/.cromwell

wb cromwell generate-config \
  --google-bucket-name=dataproc-staging-wb-cordial-diamond-9893 \
  --dir=/home/jupyter/.cromwell

cat << EOF > /home/jupyter/.cromwell/cromwell.override.conf
include "cromwell.conf"


backend.providers.GCPBATCH.config {
  concurrent-job-limit = 100
}

call-caching {
  enabled = true
}

system {
  job-rate-control {
    jobs = 10
    per = 10 seconds
  }
  workflow-heartbeats {
    ttl = 20 minutes
    write-failure-shutdown-duration = 15 minutes
    write-batch-size = 250
  }
}

engine {
  filesystems {
    gcs {
      auth = "application_default"
    }
  }
}

database {
  profile = "slick.jdbc.HsqldbProfile$"

  db {
    driver = "org.hsqldb.jdbcDriver"

    url = """
      jdbc:hsqldb:file:/home/jupyter/.cromwell/db/cromwell;
      shutdown=false;
      hsqldb.default_table_type=cached;
      hsqldb.tx=mvcc;
      hsqldb.result_max_memory_rows=10000;
      hsqldb.large_data=true;
      hsqldb.script_format=3
    """

    connectionTimeout = 120000
    numThreads = 2
  }

  insert-batch-size = 2000
}

google {
  cloud-sdk-image-url = "gcr.io/google.com/cloudsdktool/google-cloud-cli:slim"
}

backend {
  default = "GCPBATCH"
  providers {

    # Disables the Local backend
    Local.config.root = "/dev/null"

    GCPBATCH {

        actor-factory = "cromwell.backend.google.batch.GcpBatchBackendLifecycleActorFactory"

      config {
        project = "${GOOGLE_PROJECT}"
        concurrent-job-limit = 3000
        root = "gs://dataproc-temp-wb-cordial-diamond-9893/workflows/cromwell-executions"


          virtual-private-cloud {
          network-name = "projects/${GOOGLE_PROJECT}/global/networks/network"
          subnetwork-name = "projects/${GOOGLE_PROJECT}/regions/*/subnetworks/subnetwork"

        }

        batch {
          auth = "application_default"
          compute-service-account = "${PET_SA_EMAIL}"
          location = "us-central1"
        }

        default-runtime-attributes {
          noAddress: true
        }

        filesystems {
          gcs {
              auth = "application_default"
          }
        }
      }
    }
  }
}
EOF

# Ensure cromwell local database exists
DBDIR=/home/jupyter/.cromwell/db
mkdir -p $DBDIR
if [[ -d "$DBDIR" ]]; then
  size_gb=$( du -sBG "$DBDIR" | cut -f1 | tr -d 'G' )
  echo "Cromwell DB size: $size_gb GB"
fi

# Create cromshell config file
if [ ! -e /home/jupyter/.cromshell ]; then
  mkdir /home/jupyter/.cromshell
fi

cat << EOF > /home/jupyter/.cromshell/cromshell_config.json
{
  "cromwell_server": "http://localhost:8000",
  "requests_timeout": 5
}
EOF

# Launch cromwell in server mode
java \
  -Xms12G \
  -Xmx48G \
  -XX:+UseG1GC \
  -Dconfig.file=/home/jupyter/.cromwell/cromwell.override.conf \
  -Xlog:gc*:file=/home/jupyter/.cromwell/gc.log:time \
  -XX:+HeapDumpOnOutOfMemoryError \
  -XX:HeapDumpPath=/home/jupyter/.cromwell/ \
  -jar $CROMWELL_JAR \
  server
