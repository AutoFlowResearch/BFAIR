import React from 'react';
import renderer from 'react-test-renderer';
import { Provider } from 'react-redux';
import configureMockStore from 'redux-mock-store';
import Login from './login';

const mockStore = configureMockStore();
const store = mockStore({ auth: { isAuthenticating: true } });

describe('<Login />', () => {
  it('should render correctly', () => {
    const component = renderer.create(
      <Provider store={store}>
        <Login />
      </Provider>,
    );
    const tree = component.toJSON();
    expect(tree).toMatchSnapshot();
  });

  // it("should capture username onChange", () => {
  //   const component = mount(<Login />);
  //   const input = component.find("input").at(0);
  //   input.instance().value = "ketan@amrutraj.com";
  //   input.simulate("change");
  //   expect(component.state().data.username).toEqual("ketan@amrutraj.com");
  // });

  // it("should display button text loading on click", () => {
  //   const component = renderer.create(<Login />);
  //   const rootInstance = component.root;
  //   console.log(rootInstance);
  // });
});
