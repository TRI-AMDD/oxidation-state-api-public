# oxidation-state-api-public
The public repository for the Java oxidation state API

## Gradle library instructions
1. Clone the repo and add as a folder next to the existing java code.
2. In your `build.gradle` file, add as a dependency

```java
dependencies {
    implementation project(':oxidation-state-api-public')
}
```

3. Import the library as needed
```java
import tri.oxidationstates.webapi.WebOxidationAnalyzer;
    
WebOxidationAnalyzer analyzer = new WebOxidationAnalyzer(paramFileName, polyIonDir);
```