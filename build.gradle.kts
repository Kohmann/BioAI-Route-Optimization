import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

plugins {
    //kotlin("jvm") version "1.5.10"
    application
    kotlin("jvm") version "1.6.10" // or kotlin("multiplatform") or any other kotlin plugin
    kotlin("plugin.serialization") version "1.6.10"
}

group = "me.kohmis"
version = "1.0-SNAPSHOT"

repositories {
    mavenCentral()
}

dependencies {
    testImplementation(kotlin("test"))
    implementation("org.jetbrains.kotlinx:kotlinx-serialization-json:1.3.2")
    implementation("com.github.haifengl:smile-kotlin:2.6.0")
    implementation("org.slf4j:slf4j-simple:1.7.36")
}

tasks.test {
    useJUnitPlatform()
}

tasks.withType<KotlinCompile> {
    kotlinOptions.jvmTarget = "1.8"
}

application {
    mainClass.set("MainKt")
}