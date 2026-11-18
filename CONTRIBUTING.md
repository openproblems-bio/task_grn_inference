# Contribution guidelines

We encourage contributions from the community. To contribute:

* Star this repository: Click the star button in the top-right corner of the repository to show your support.
* Fork the repository: Start by forking this repository to your account.
* Develop your component: Create your Viash component, ensuring it aligns with our best practices (detailed below).
* Submit a pull request: After testing your component, submit a pull request for review.

## Installation

You need to have Docker, Java, and Viash installed. Follow
[these instructions](https://openproblems.bio/documentation/fundamentals/requirements)
to install the required dependencies.

## Getting started

### Cloning the repository

To get started, fork the repository and clone it to your local machine:

```bash
git clone <git url to fork> --recursive

cd <name of the repository>
```

NOTE: If you forgot to clone the repository with the `--recursive` flag, you can run the following command to update the submodules:

```bash
git submodule update --init --recursive
```

### Downloading the test resources

Next, you should download the test resources:

```bash
scripts/sync_resources.sh
```

You may need to run this script again if the resources are updated.

## Good first steps

### Creating a new method

To create a new method, you can use the following command:

```bash
# in Python:
common/scripts/create_component \
  --name my_method \
  --language python \
  --type method

# or in R:
common/scripts/create_component \
  --name my_method \
  --language r \
  --type method
```

This will create a new method in `src/methods/my_method`. Next, you should:

* Fill in the component's metadata
* Specify dependencies
* Implement the method's code
* Run the unit test

Please review our documentation on [creating a new method](https://openproblems.bio/documentation/create_component/add-a-method/) for more information on how to do this.


### Creating a new metric

Creating a new metric is similar to creating a new method. You can use the following command:

```bash
# in Python:
common/scripts/create_component \
  --name my_metric \
  --language python \
  --type metric

# or in R:
common/scripts/create_component \
  --name my_metric \
  --language r \
  --type metric
```

This will create a new metric in `src/metrics/my_metric`. Next, you should:

* Fill in the component's metadata
* Specify dependencies
* Implement the metric's code
* Run the unit test

Please review our documentation on [creating a new metric](https://openproblems.bio/documentation/create_component/add-a-metric/) for more information.


## Frequently used commands

To get started, you can run the following commands:

### Sync the test data

To sync the test data, run the following command:

```bash
scripts/sync_resources.sh
```

### Building the components

To run the benchmark, you first need to build the components.

```bash
viash ns build --parallel --setup cachedbuild
```

For each of the components, this will:

* Build a Docker image
* Build an executable at `target/executable/<component_name>`
* Build a Nextflow module at `target/nextflow/<component_name>`

### Running the unit tests

To run the unit test of one component, you can use the following command:

```bash
viash test src/path/to/config.vsh.yaml
```

To run the unit tests for all components, you can use the following command:

```bash
viash ns test --parallel
```

### Running the benchmark

To run the benchmark on the test resources, you can use the following command:

```bash
scripts/run_benchmark/run_test_local.sh
```
