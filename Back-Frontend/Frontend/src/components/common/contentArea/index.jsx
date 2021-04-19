import React, { Component } from "react";
import { Redirect, Route, Switch } from "react-router-dom";
// import { Route } from "../protectedRoutes";
import Dashboard from "../../pages/admin/dashboard";
import UserForm from "../../pages/admin/NewUser";
import Posts from "../../pages/admin/posts";
import NewPost from "../../pages/admin/posts/NewPost";
import NewRecipe from "../../pages/admin/recipe/NewRecipe";
import RecipeList from "../../pages/admin/recipe/RecipeList";
import SimpleForm from "../../pages/admin/SimpleForm";
import SampleSubmission from "../../pages/admin/userManagement"

class ContentArea extends Component {
  state = {
    role: "admin"
  };

  componentDidMount() {
    // console.log(this.props.location);
    // this.props.history.push(`/mainpage/${this.state.role}/dashboard`);
  }

  render() {
    return (
      <div className="content-wrapper">
        <Switch>
          <Route
            path="/mainpage/admin/dashboard"
            roles={["admin"]}
            component={Dashboard}
          />

          <Route
            path="/mainpage/admin/user-management/:id"
            component={UserForm}
          />

          <Route
            path="/mainpage/admin/doe-module"
            component={SampleSubmission}
          />

          <Route
            path="/mainpage/admin/posts/:id"
            component={NewPost}
          />

          <Route path="/mainpage/admin/posts" component={Posts} />

          <Route
            path="/mainpage/admin/simple-form"
            component={SimpleForm}
          />
          <Route
            path="/mainpage/admin/recipe/:id"
            component={NewRecipe}
          />
          <Route
            path="/mainpage/admin/recipe"
            component={RecipeList}
          />

          <Redirect to="/mainpage/admin/dashboard" />
        </Switch>
      </div>
    );
  }
}

export default ContentArea;
