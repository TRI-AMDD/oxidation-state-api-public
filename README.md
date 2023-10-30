# oxidation-state-api-public
The public repository for the Java library used by the [oxidation state analyzer](https://oxi.matr.io).

## How to build the Java library
1. Install [gradle](https://gradle.org/install/).
2. Clone the repo.
3. From the command line, run:

```shell
gradle clean build javadoc
```
- The library will be in build/lib.
- The dependencies will be in build/dependencies.
- The javadoc will be in build/docs.

## How to call the library from Python
1. Build the library as described above.
2. Install [jPype](https://jpype.readthedocs.io/).
3. Follow the example in the examples/python/analyze_oxidation_states.ipynb notebook.

## How to include the library in a Java project using Gradle
1. Clone the repo and add it as a folder next to the existing java code.
2. In your `build.gradle` file, add as the oxidation state API as a dependency:
        
```java
dependencies {
    implementation project(':oxidation-state-api-public')
}
```
       